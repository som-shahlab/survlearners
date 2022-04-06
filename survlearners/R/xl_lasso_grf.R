#' @include utils.R
#'
#' @title X-learner of lasso and grf
#'
#' @description  X-learner, implemented via grf (for nuisance parameter estimation)
#' and glmnet (lasso) with 'coxph' distribution (for target function CATE estimation)
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
#' @param ps The propensity score
#' @param cen_fit The choice of model fitting for censoring
#' @examples
#' \dontrun{
#' n = 1000; p = 25
#' times = 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[,1] + (-0.5 - 1 * X[,2]) * W))^2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC/(4^2))^(1/2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' data <- list(X = X, W = W, Y = Y, D = D)
#' data.test <- list(X = X, W = W, Y = Y, D = D)
#'
#' xlas_grf_surv_cate = estimate_ipcw_las_grf_xl(data, data.test, times, ps = 0.5)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
estimate_ipcw_las_grf_xl <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){
Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
times.index <- findInterval(times, Y.grid)

# fit model on W==1
grffit1 <- grf::survival_forest(data$X[data$W==1,],
                                data$Y[data$W==1],
                                data$D[data$W==1],
                                alpha = alpha,
                                prediction.type = "Nelson-Aalen",
                                failure.times = Y.grid)
surf1 <- rep(NA, length(data$W))
surf1[data$W==1] <- predict(grffit1)$predictions[, times.index]
surf1[data$W==0] <- predict(grffit1, data$X[data$W==0,])$predictions[, times.index]

# fit model on W==0
grffit0 <- grf::survival_forest(data$X[data$W==0,],
                                data$Y[data$W==0],
                                data$D[data$W==0],
                                alpha = alpha,
                                prediction.type = "Nelson-Aalen",
                                failure.times = Y.grid)
surf0 <- rep(NA, length(data$W))
surf0[traindat$W==0] <- predict(grffit0)$predictions[, times.index]
surf0[traindat$W==1] <- predict(grffit0, data$X[data$W==1,])$predictions[, times.index]

Tgrf1 <- 1-surf1
Tgrf0 <- 1-surf0

# IPCW weights
if(cen_fit == "KM"){
  shuffle <- sample(nrow(traindat))
  kmdat <- data.frame(Y = data$Y[shuffle], D = data$D[shuffle])
  folds <- cut(seq(1, nrow(kmdat)), breaks=10, labels=FALSE)
  C.Y.hat <- rep(NA, nrow(kmdat))
  for(z in 1:10){
    testIndexes <- which(folds==z, arr.ind=TRUE)
    testData <- kmdat[testIndexes, ]
    trainData <- kmdat[-testIndexes, ]
    km_fit <- survfit(Surv(trainData$Y, 1 - trainData$D) ~ 1)
    cent <- testData$Y
    cent[testData$D==0] <- times
    C.Y.hat[testIndexes] <- summary(km_fit, times = cent)$surv
  }
  shudat <- data.frame(shuffle, C.Y.hat)
  C.Y.hat <- shudat[order(shuffle), ]$C.Y.hat
}else if (cen_fit == "survival.forest"){
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  c_fit <- grf::survival_forest(cbind(data$W, data$X),
                                data$Y,
                                1 - data$D,
                                failure.times = Y.grid,
                                alpha = alpha,
                                prediction.type = "Nelson-Aalen")
  C.hat <- predict(c_fit)$predictions
  cent <- data$Y; cent[data$D==0] <- times
  cen.times.index <- findInterval(cent, c_fit$failure.times)
  C.Y.hat <- C.hat[cbind(1:length(data$Y), cen.times.index)]
}
ipcw <- 1 / C.Y.hat

# Propensity score
if (is.null(ps)){
  stop("propensity score needs to be supplied")
}else{
  ps.train <- rep(ps, length(data$W))
  ps_test <- rep(ps, length(data.test$W))
}

weight <- ipcw/ps.train  # censoring weight * treatment weight

# X-learner
tempdat <- data.frame(Y = data$Y, D = data$D, W = data$W, weight, X = data$X, Tgrf0, Tgrf1)
binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]
binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0
binary_data <- binary_data[complete.cases(binary_data), ]
b_data <- list(Y = binary_data$Y, D = binary_data$D, W = binary_data$W, X = binary_data$X,
               wt = binary_data$weight, mu0 = binary_data$Tgrf0, mu1 = binary_data$Tgrf1)

foldid <- sample(rep(seq(10), length = length(b_data$Y[b_data$W==1])))
XLfit1 <- glmnet::cv.glmnet(b_data$X[b_data$W==1, ],
                            b_data$D[b_data$W==1] - b_data$mu0[b_data$W==1],
                            weights = b_data$wt[b_data$W==1],
                            foldid = foldid,
                            alpha = 1)
XLtau1 <- as.vector(-predict(XLfit1, data.test$X, s = "lambda.min"))

foldid <- sample(rep(seq(10), length = length(b_data$Y[b_data$W==0])))
XLfit0 <- glmnet::cv.glmnet(b_data$X[b_data$W==0, ],
                            b_data$mu1[b_data$W==0] - b_data$D[b_data$W==0],
                            weights = b_data$wt[b_data$W==0],
                            foldid = foldid,
                            alpha = 1)
XLtau0 <- as.vector(-predict(XLfit0, data.test$X, s = "lambda.min"))

# weighted CATE
pred_X_las_grf <- XLtau1 * (1 - ps_test) + XLtau0 * ps_test
as.vector(pred_X_las_grf)
}
