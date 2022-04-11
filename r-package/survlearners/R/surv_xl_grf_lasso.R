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
#' @param cen.fit The choice of model fitting for censoring
#' @param newX The test data set (covariates only)
#' @examples
#' \donttest{
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
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.xl.grf.lasso.fit = surv_xl_grf_lasso(X, W, Y, D, times, ps = 0.5)
#' cate = predict(surv.xl.grf.lasso.fit)
#' cate.test = predict(surv.xl.grf.lasso.fit, X.test)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_xl_grf_lasso <- function(X, W, Y, D, times, alpha = 0.05, ps = NULL, cen.fit = "KM"){

# fit model on W==1
grffit1 <- grf::survival_forest(X[W==1,],
                                Y[W==1],
                                D[W==1],
                                alpha = alpha,
                                prediction.type = "Nelson-Aalen")
surf1 <- rep(NA, length(W))
times.index <- findInterval(times, grffit1$failure.times)
surf1[W==1] <- predict(grffit1)$predictions[, times.index]
surf1[W==0] <- predict(grffit1, X[W==0,])$predictions[, times.index]

# fit model on W==0
grffit0 <- grf::survival_forest(X[W==0,],
                                Y[W==0],
                                D[W==0],
                                alpha = alpha,
                                prediction.type = "Nelson-Aalen")
surf0 <- rep(NA, length(W))
times.index <- findInterval(times, grffit0$failure.times)
surf0[W==0] <- predict(grffit0)$predictions[, times.index]
surf0[W==1] <- predict(grffit0, X[W==1,])$predictions[, times.index]

Tgrf1 <- 1-surf1
Tgrf0 <- 1-surf0

# IPCW weights
if(cen.fit == "KM"){
  shuffle <- sample(length(Y))
  kmdat <- data.frame(Y = Y[shuffle], D = D[shuffle])
  folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
  C.Y.hat <- rep(NA, nrow(kmdat))
  for(z in 1:10){
    testIndexes <- which(folds==z, arr.ind=TRUE)
    testData <- kmdat[testIndexes, ]
    trainData <- kmdat[-testIndexes, ]
    c.fit <- survival::survfit(survival::Surv(trainData$Y, 1 - trainData$D) ~ 1)
    cent <- testData$Y; cent[testData$D==0] <- times
    C.Y.hat[testIndexes] <- summary(c.fit, times = cent)$surv
  }
  shudat <- data.frame(shuffle, C.Y.hat)
  C.Y.hat <- shudat[order(shuffle), ]$C.Y.hat
}else if (cen.fit == "survival.forest"){
  c.fit <- grf::survival_forest(cbind(W, X),
                                Y,
                                1 - D,
                                alpha = alpha,
                                prediction.type = "Nelson-Aalen")
  C.hat <- predict(c.fit)$predictions
  cent <- Y; cent[D==0] <- times
  cen.times.index <- findInterval(cent, c.fit$failure.times)
  C.Y.hat <- C.hat[cbind(1:length(Y), cen.times.index)]
}
ipcw <- 1 / C.Y.hat

# Propensity score
if (is.null(ps)){
  stop("propensity score needs to be supplied")
}else{
  ps.train <- rep(ps, length(W))
}

weight <- ipcw/ps.train  # censoring weight * treatment weight

# X-learner
tempdat <- data.frame(Y = Y, D = D, W = W, weight, X, Tgrf0, Tgrf1)
binary.data <- tempdat[tempdat$D==1|tempdat$Y > times,]
binary.data$D[binary.data$D==1 & binary.data$Y > times] <- 0
binary.data <- binary.data[complete.cases(binary.data), ]
b.data <- list(Y = binary.data$Y, D = binary.data$D, W = binary.data$W,
               X = as.matrix(binary.data[,5:(ncol(binary.data)-2)]),
               wt = binary.data$weight, mu0 = binary.data$Tgrf0, mu1 = binary.data$Tgrf1)

foldid <- sample(rep(seq(10), length = length(b.data$Y[b.data$W==1])))
XLfit1 <- glmnet::cv.glmnet(b.data$X[b.data$W==1, ],
                            b.data$D[b.data$W==1] - b.data$mu0[b.data$W==1],
                            weights = b.data$wt[b.data$W==1],
                            foldid = foldid,
                            alpha = 1)
XLtau1 <- as.vector(-predict(XLfit1, X, s = "lambda.min"))

foldid <- sample(rep(seq(10), length = length(b.data$Y[b.data$W==0])))
XLfit0 <- glmnet::cv.glmnet(b.data$X[b.data$W==0, ],
                            b.data$mu1[b.data$W==0] - b.data$D[b.data$W==0],
                            weights = b.data$wt[b.data$W==0],
                            foldid = foldid,
                            alpha = 1)
XLtau0 <- as.vector(-predict(XLfit0, X, s = "lambda.min"))

# weighted CATE
pred.X.las.grf <- as.vector(XLtau1 * (1 - ps.train) + XLtau0 * ps.train)

ret <- list(fit1 = XLfit1,
            fit0 = XLfit0,
            ps = ps.train,
            tau = pred.X.las.grf)
class(ret) <- 'surv_xl_grf_lasso'
ret
}

#' predict for surv_xl_grf_lasso
#'
#' get estimated tau(X) using the trained surv_xl_grf_lasso model
#'
#' @param object An surv_xl_grf_lasso object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param ps The propensity score
#' @param ... Additional arguments (currently not used)
#'
#' @examples
#' \donttest{
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
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.xl.grf.lasso.fit = surv_xl_grf_lasso(X, W, Y, D, times, ps = 0.5)
#' cate = predict(surv.xl.grf.lasso.fit)
#' cate.test = predict(surv.xl.grf.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_xl_grf_lasso = function(object,
                                     newdata = NULL,
                                     ps = NULL,
                                     ...) {
  if(is.null(newdata)){
    return(object$tau)
  }else{
    XLtau1 <- as.vector(-predict(object$fit1, newdata, s = "lambda.min"))
    XLtau0 <- as.vector(-predict(object$fit0, newdata, s = "lambda.min"))
    if(is.null(ps)){
      ps <- rep(object$ps[1], nrow(newdata))
    }else{
      if(length(ps)==1){
      ps <- rep(ps, nrow(newdata))
      }
    }
    return(as.vector(XLtau1 * (1 - ps) + XLtau0 * ps))
  }
}
