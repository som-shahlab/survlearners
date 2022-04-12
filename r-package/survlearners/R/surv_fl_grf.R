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
#' @param W.hat The propensity score
#' @param cen.fit The choice of model fitting for censoring
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
#' surv.fl.grf.fit = surv_fl_grf(X, Y, W, D, times, W.hat = 0.5)
#' cate = predict(surv.fl.grf.fit)
#' cate.test = predict(surv.fl.grf.fit, X.test)
#' }
#' @return A surv_fl_grf object
#' @export
surv_fl_grf <- function(X, Y, W, D, times, alpha = 0.05, W.hat = NULL, cen.fit = "Kaplan-Meier"){

  # IPCW weights
  if (cen.fit == "Kaplan-Meier"){
    shuffle <- sample(length(Y))
    kmdat <- data.frame(Y = Y[shuffle], D = D[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
    C.hat <- rep(NA, nrow(kmdat))
    for(z in 1:10){
      testIndexes <- which(folds==z, arr.ind=TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      c.fit <- survival::survfit(survival::Surv(trainData$Y, 1 - trainData$D) ~ 1)
      cent <- testData$Y; cent[testData$D==0] <- times
      C.hat[testIndexes] <- summary(c.fit, times = cent)$surv
    }
    shudat <- data.frame(shuffle, C.hat)
    C.hat <- shudat[order(shuffle), ]$C.hat
  } else if (cen.fit == "survival.forest"){
    c.fit <- grf::survival_forest(cbind(W, X),
                                  Y,
                                  1 - D,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
    C.hat <- predict(c.fit)$predictions
    cent <- Y; cent[D==0] <- times
    cen.times.index <- findInterval(cent, c.fit$failure.times)
    C.hat <- C.hat[cbind(1:length(Y), cen.times.index)]
  }
  sample.weights <- 1 / C.hat

  # Propensity score
  if (is.null(W.hat)){
    stop("propensity score needs to be supplied")
  } else {
    W.hat <- rep(W.hat, length(Y))
  }

  # Subset of uncensored subjects
  tempdat <- data.frame(Y = Y, D = D, W = W, W.hat, sample.weights, X)
  binary.data <- tempdat[tempdat$D==1|tempdat$Y > times,]
  binary.data$D[binary.data$D==1 & binary.data$Y > times] <- 0
  binary.data <- binary.data[complete.cases(binary.data), ]
  b.data <- list(Y = binary.data$Y, D = binary.data$D, W = binary.data$W,
                 X = as.matrix(binary.data[,6:ncol(binary.data)]),
                 sample.weights = binary.data$sample.weights, W.hat = binary.data$W.hat)

  Z <- b.data$W * b.data$D / b.data$W.hat - (1 - b.data$W) * b.data$D / (1 - b.data$W.hat)
  tau.fit <- grf::regression_forest(b.data$X, Z, sample.weights = b.data$sample.weights)
  tau.hat <- -predict(tau.fit, X)

  ret <- list(tau.fit = tau.fit,
              tau.hat = tau.hat)
  class(ret) <- "surv_fl_grf"
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
#' surv.fl.grf.fit = surv_fl_grf(X, Y, W, D, times, W.hat = 0.5)
#' cate = predict(surv.fl.grf.fit)
#' cate.test = predict(surv.fl.grf.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_fl_grf = function(object,
                               newdata = NULL,
                               ...) {
  if (is.null(newdata)){
    return(object$tau.hat)
  } else {
    return(-predict(object$tau.fit, newdata))
  }
}
