#' @title S-learner of grf
#'
#' @description  S-learner, implemented via survival_forest (grf package)
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' times <- 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[ ,1] + (-0.5 - 1 * X[ ,2]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.sl.grf.fit <- surv_sl_grf(X, Y, W, D, times)
#' cate <- predict(surv.sl.grf.fit)
#' cate.test <- predict(surv.sl.grf.fit, X.test)
#' }
#' @return A surv_sl_grf object
#' @export
surv_sl_grf <- function(X, Y, W, D, times, alpha = 0.05) {

  tau.fit <- grf::survival_forest(cbind(W, X),
                                  Y,
                                  D,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")

  index <- findInterval(times, tau.fit$failure.times)
  surf1 <- predict(tau.fit, cbind(rep(1, nrow(X)), X))$predictions[ ,index]
  surf0 <- predict(tau.fit, cbind(rep(0, nrow(X)), X))$predictions[ ,index]
  tau.hat <- surf1 - surf0

  ret <- list(tau.fit = tau.fit,
              tau.hat = tau.hat,
              times = times)
  class(ret) <- "surv_sl_grf"
  ret
}

#' predict for surv_sl_grf
#'
#' get estimated tau(X) using the trained surv_sl_grf model
#'
#' @param object An surv_sl_grf object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param times The prediction time of interest
#' @param ... Additional arguments (currently not used)
#'
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' times <- 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[ ,1] + (-0.5 - 1 * X[ ,2]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.sl.grf.fit <- surv_sl_grf(X, Y, W, D, times)
#' cate <- predict(surv.sl.grf.fit)
#' cate.test <- predict(surv.sl.grf.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_sl_grf <- function(object,
                                newdata = NULL,
                                times = NULL,
                                ...) {
  if (is.null(newdata)) {
    return(object$tau.hat)
  } else {
    if (is.null(times)) {
      index <- findInterval(object$times, object$fit$failure.times)
    } else {
      index <- findInterval(times, object$fit$failure.times)
    }
    surf1 <- predict(object$tau.fit, cbind(rep(1, nrow(newdata)), newdata))$predictions[ ,index]
    surf0 <- predict(object$tau.fit, cbind(rep(0, nrow(newdata)), newdata))$predictions[ ,index]
    return(surf1 - surf0)
  }
}
