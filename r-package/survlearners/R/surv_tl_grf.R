#' @title T-learner of grf
#'
#' @description  T-learner, implemented via survival_forest (grf package)
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' t0 <- 0.2
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
#' surv.tl.grf.fit <- surv_tl_grf(X, Y, W, D, t0)
#' cate <- predict(surv.tl.grf.fit)
#' cate.test <- predict(surv.tl.grf.fit, X.test)
#' }
#' @return A surv_tl_grf object
#' @export
surv_tl_grf <- function(X, Y, W, D, t0, alpha = 0.05) {
  # Model for W = 1
  grffit1 <- grf::survival_forest(X[W == 1, ],
                                  Y[W == 1],
                                  D[W == 1],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
  index <- findInterval(t0, grffit1$failure.times)
  if (index == 0) {
    surf1 <- rep(1, nrow(X))
  } else {
    surf1 <- predict(grffit1, X)$predictions[ ,index]
  }

  # Model for W = 0
  grffit0 <- grf::survival_forest(X[W == 0, ],
                                  Y[W == 0],
                                  D[W == 0],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
  index <- findInterval(t0, grffit0$failure.times)
  if (index == 0) {
    surf0 <- rep(1, nrow(X))
  } else {
    surf0 <- predict(grffit0, X)$predictions[ ,index]
  }

  tau.hat <- surf1 - surf0

  ret <- list(fit1 = grffit1,
              fit0 = grffit0,
              tau.hat = tau.hat,
              t0 = t0)
  class(ret) <- "surv_tl_grf"
  ret
}


#' predict for surv_tl_grf
#'
#' get estimated tau(X) using the trained surv_tl_grf model
#'
#' @param object An surv_tl_grf object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param t0 The prediction time of interest
#' @param ... Additional arguments (currently not used)
#'
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' t0 <- 0.2
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
#' surv.tl.grf.fit <- surv_tl_grf(X, Y, W, D, t0)
#' cate <- predict(surv.tl.grf.fit)
#' cate.test <- predict(surv.tl.grf.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_tl_grf <- function(object,
                                newdata = NULL,
                                t0 = NULL,
                                ...) {
  if (is.null(newdata)) {
    return(object$tau.hat)
  } else {
    if (is.null(t0)) {
      index1 <- findInterval(object$t0, object$fit1$failure.times)
      index0 <- findInterval(object$t0, object$fit0$failure.times)
    } else {
      index1 <- findInterval(t0, object$fit1$failure.times)
      index0 <- findInterval(t0, object$fit0$failure.times)
    }
    if (index1 == 0) {
      surf1 <- rep(1, nrow(newdata))
    } else {
      surf1 <- predict(object$fit1, newdata)$predictions[ ,index1]
    }

    if (index0 == 0) {
      surf0 <- rep(1, nrow(newdata))
    } else {
      surf0 <- predict(object$fit0, newdata)$predictions[ ,index0]
    }

    return(surf1 - surf0)
  }
}
