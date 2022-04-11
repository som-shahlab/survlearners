#' @title F-learner of grf
#'
#' @description  F-learner, implemented via survival_forest in the grf package
#'
#' @param X The baseline covariates
#' @param W The treatment variable (0 or 1)
#' @param Y The follow-up time
#' @param D The event indicator
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
#' @param ps The propensity score
#' @param cen_fit The choice of model fitting for censoring
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
#' surv_fl_grf_fit = surv_fl_grf(X, W, Y, D, times, ps = 0.5)
#' cate = predict(surv_fl_grf_fit)
#' cate.test = predict(surv_fl_grf_fit, X.test)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_fl_grf <- function(X, W, Y, D, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){

  # IPCW weights
  if(cen_fit == "KM"){
    shuffle <- sample(length(Y))
    kmdat <- data.frame(Y = Y[shuffle], D = D[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
    C.Y.hat <- rep(NA, nrow(kmdat))
    for(z in 1:10){
      testIndexes <- which(folds==z, arr.ind=TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      c_fit <- survival::survfit(survival::Surv(trainData$Y, 1 - trainData$D) ~ 1)
      cent <- testData$Y; cent[testData$D==0] <- times
      C.Y.hat[testIndexes] <- summary(c_fit, times = cent)$surv
    }
    shudat <- data.frame(shuffle, C.Y.hat)
    C.Y.hat <- shudat[order(shuffle), ]$C.Y.hat
  }else if (cen_fit == "survival.forest"){
    c_fit <- grf::survival_forest(cbind(W, X),
                                  Y,
                                  1 - D,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
    C.hat <- predict(c_fit)$predictions
    cent <- Y; cent[D==0] <- times
    cen.times.index <- findInterval(cent, c_fit$failure.times)
    C.Y.hat <- C.hat[cbind(1:length(Y), cen.times.index)]
  }
  ipcw <- 1 / C.Y.hat

  # Propensity score
  if (is.null(ps)){
    stop("propensity score needs to be supplied")
  }else{
    pscore <- rep(ps, length(Y))
  }

  # Subset of uncensored subjects
  tempdat <- data.frame(Y = Y, D = D, W = W, pscore, ipcw, X)
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0
  binary_data <- binary_data[complete.cases(binary_data), ]
  b_data <- list(Y = binary_data$Y, D = binary_data$D, W = binary_data$W,
                 X = as.matrix(binary_data[,6:ncol(binary_data)]),
                 wt = binary_data$ipcw, ps = binary_data$pscore)

  Z <- b_data$W * b_data$D / b_data$ps - (1 - b_data$W) * b_data$D / (1 - b_data$ps)
  fgrf_fit <- grf::regression_forest(b_data$X, Z, sample.weights = b_data$wt)
  fgrf_tau <- -predict(fgrf_fit, X)

  ret <- list(fit = fgrf_fit,
              tau = fgrf_tau)
  class(ret) <- 'surv_fl_grf'
  ret
}

#' predict for surv_fl_grf
#'
#' get estimated tau(X) using the trained surv_fl_grf model
#'
#' @param object An surv_fl_grf object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
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
#' surv_fl_grf_fit = surv_fl_grf(X, W, Y, D, times, ps = 0.5)
#' cate = predict(surv_fl_grf_fit)
#' cate.test = predict(surv_fl_grf_fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_fl_grf = function(object,
                               newdata = NULL,
                               ...) {
  if(is.null(newdata)){
    return(object$tau)
  }else{
    return(-predict(object$fit, newdata))
  }
}
