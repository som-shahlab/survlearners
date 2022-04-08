#' @title F-learner of lasso
#'
#' @description  F-learner, implemented via glmnet (lasso) with 'coxph' distribution
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
#' flasso_surv_cate = estimate_ipcw_lasso_fl(data, data.test, times, ps = 0.5)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
estimate_ipcw_lasso_fl <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){

  # IPCW weights
  if(cen_fit == "KM"){
    shuffle <- sample(length(data$Y))
    kmdat <- data.frame(Y = data$Y[shuffle], D = data$D[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
    c_hat <- rep(NA, nrow(kmdat))
    for(z in 1:10){
      testIndexes <- which(folds==z, arr.ind=TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      c_fit <- survfit(Surv(trainData$Y, 1 - trainData$D) ~ 1)
      cent <- testData$Y
      cent[testData$D==0] <- times
      c_hat[testIndexes] <- summary(c_fit, times = cent)$surv
    }
    shudat <- data.frame(shuffle, c_hat)
    c_hat <- shudat[order(shuffle), ]$c_hat
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
    c_hat <- C.hat[cbind(1:length(data$Y), cen.times.index)]
  }
  ipcw <- 1 / c_hat

  # Propensity score
  if (is.null(ps)){
    stop("propensity score needs to be supplied")
  }else{
    ps_score <- rep(ps, length(data$Y))
  }

  # Subset of uncensored subjects
  tempdat <- data.frame(Y = data$Y, D = data$D, W = data$W, ps_score, ipcw, X = data$X)
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0
  binary_data <- binary_data[complete.cases(binary_data), ]
  b_data <- list(Y = binary_data$Y, D = binary_data$D, W = binary_data$W, X = binary_data$X,
                 wt = binary_data$ipcw, ps = binary_data$ps_score)

  flasso_fit <- Flasso(x = b_data$X,
                       tx = b_data$W,
                       y = b_data$D,
                       pscore = b_data$ps,
                       weight = b_data$wt)
  pred_flasso <- as.vector(-predict(flasso_fit, data.test$X))
  pred_flasso
}
