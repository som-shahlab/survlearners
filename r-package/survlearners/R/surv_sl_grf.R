#' @title S-learner of grf
#'
#' @description  S-learner, implemented via survival_forest (grf package)
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param args.grf.nuisance Input arguments for a grf model that estimates nuisance parameters
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' t0 <- 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[ ,1, drop = FALSE] + (-0.5 - 1 * X[ ,2, drop = FALSE]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.sl.grf.fit <- surv_sl_grf(X, Y, W, D, t0)
#' cate <- predict(surv.sl.grf.fit)
#' cate.test <- predict(surv.sl.grf.fit, X.test)
#' }
#' @return A surv_sl_grf object
#' @export
surv_sl_grf <- function(X, Y, W, D, t0, args.grf.nuisance = list()) {

  args.grf.nuisance <- list(failure.times = NULL,
                            num.trees = max(50, 2000 / 4),
                            sample.weights = NULL,
                            clusters = NULL,
                            equalize.cluster.weights = FALSE,
                            sample.fraction = 0.5,
                            mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                            min.node.size = 15,
                            honesty = TRUE,
                            honesty.fraction = 0.5,
                            honesty.prune.leaves = TRUE,
                            alpha = 0.05,
                            prediction.type = "Nelson-Aalen",
                            compute.oob.predictions = TRUE,
                            num.threads = NULL,
                            seed = runif(1, 0, .Machine$integer.max))

  tau.fit <- do.call(grf::survival_forest, c(list(X = cbind(W, X), Y = Y, D = D), args.grf.nuisance))
  index <- findInterval(t0, tau.fit$failure.times)
  if (index == 0) {
    surf1 <- rep(1, nrow(X))
    surf0 <- rep(1, nrow(X))
  } else {
    surf1 <- predict(tau.fit, cbind(rep(1, nrow(X)), X))$predictions[ ,index, drop = FALSE]
    surf0 <- predict(tau.fit, cbind(rep(0, nrow(X)), X))$predictions[ ,index, drop = FALSE]
  }

  tau.hat <- surf1 - surf0

  ret <- list(tau.fit = tau.fit,
              tau.hat = tau.hat,
              t0 = t0)
  class(ret) <- "surv_sl_grf"
  ret
}

#' predict for surv_sl_grf
#'
#' get estimated tau(X) using the trained surv_sl_grf model
#'
#' @param object An surv_sl_grf object
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
#' T <- (numeratorT / exp(1 * X[ ,1, drop = FALSE] + (-0.5 - 1 * X[ ,2, drop = FALSE]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.sl.grf.fit <- surv_sl_grf(X, Y, W, D, t0)
#' cate <- predict(surv.sl.grf.fit)
#' cate.test <- predict(surv.sl.grf.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_sl_grf <- function(object,
                                newdata = NULL,
                                t0 = NULL,
                                ...) {
  if (is.null(newdata)) {
    return(object$tau.hat)
  } else {
    if (is.null(t0)) {
      index <- findInterval(object$t0, object$tau.fit$failure.times)
    } else {
      index <- findInterval(t0, object$tau.fit$failure.times)
    }

    if (index == 0) {
      surf1 <- rep(1, nrow(newdata))
      surf0 <- rep(1, nrow(newdata))
    } else {
      surf1 <- predict(object$tau.fit, cbind(rep(1, nrow(newdata)), newdata))$predictions[ ,index, drop = FALSE]
      surf0 <- predict(object$tau.fit, cbind(rep(0, nrow(newdata)), newdata))$predictions[ ,index, drop = FALSE]
    }

    return(surf1 - surf0)
  }
}
