#' @title T-learner of grf
#'
#' @description  T-learner, implemented via survival_forest (grf package)
#'
#' @param X The baseline covariates
#' @param W The treatment variable (0 or 1)
#' @param Y The follow-up time
#' @param D The event indicator
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
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
#' surv.tl.grf.fit = surv_tl_grf(X, Y, W, D, times)
#' cate = predict(surv.tl.grf.fit)
#' cate.test = predict(surv.tl.grf.fit, X.test)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_tl_grf <- function(X, Y, W, D, times, alpha = 0.05){
  # Model for W = 1
  grffit1 <- grf::survival_forest(X[W==1,],
                                  Y[W==1],
                                  D[W==1],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
  index <- findInterval(times, grffit1$failure.times)
  surf1 <- predict(grffit1, X)$predictions[, index]

  # Model for W = 0
  grffit0 <- grf::survival_forest(X[W==0,],
                                  Y[W==0],
                                  D[W==0],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
  index <- findInterval(times, grffit0$failure.times)
  surf0 <- predict(grffit0, X)$predictions[, index]

  tau.hat <- surf1 - surf0

  ret <- list(fit1 = grffit1,
              fit0 = grffit0,
              tau.hat = tau.hat,
              times = times)
  class(ret) <- 'surv_tl_grf'
  ret
}


#' predict for surv_tl_grf
#'
#' get estimated tau(X) using the trained surv_tl_grf model
#'
#' @param object An surv_tl_grf object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param times The prediction time of interest
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
#' surv.tl.grf.fit = surv_tl_grf(X, Y, W, D, times)
#' cate = predict(surv.tl.grf.fit)
#' cate.test = predict(surv.tl.grf.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_tl_grf = function(object,
                               newdata = NULL,
                               times = NULL,
                               ...) {
  if(is.null(newdata)){
    return(object$tau.hat)
  }else{
    if(is.null(times)){
      index1 <- findInterval(object$times, object$fit1$failure.times)
      index0 <- findInterval(object$times, object$fit0$failure.times)
    }else{
      index1 <- findInterval(times, object$fit1$failure.times)
      index0 <- findInterval(times, object$fit0$failure.times)
    }
    surf1 <- predict(object$fit1, newdata)$predictions[, index1]
    surf0 <- predict(object$fit0, newdata)$predictions[, index0]
    return(surf1 - surf0)
  }
}
