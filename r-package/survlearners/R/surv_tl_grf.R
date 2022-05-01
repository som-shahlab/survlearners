#' @title T-learner of grf
#'
#' @description  T-learner, implemented via survival_forest (grf package)
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param new.args.grf.nuisance Input arguments for a grf model that estimates nuisance parameters
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
#' surv.tl.grf.fit <- surv_tl_grf(X, Y, W, D, t0)
#' cate <- predict(surv.tl.grf.fit)
#' cate.test <- predict(surv.tl.grf.fit, X.test)
#' }
#' @return A surv_tl_grf object
#' @export
surv_tl_grf <- function(X, Y, W, D, t0, new.args.grf.nuisance = list()) {

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
  args.grf.nuisance[names(new.args.grf.nuisance)] <- new.args.grf.nuisance

  # Model for W = 1
  grffit1 <- do.call(grf::survival_forest, c(list(X = X[W == 1,, drop = FALSE], Y = Y[W == 1], D = D[W == 1]), args.grf.nuisance))
  surf1 <- predict(grffit1, X, failure.times = t0)$predictions

  # Model for W = 0
  grffit0 <- do.call(grf::survival_forest, c(list(X = X[W == 0,, drop = FALSE], Y = Y[W == 0], D = D[W == 0]), args.grf.nuisance))
  surf0 <- predict(grffit0, X, failure.times = t0)$predictions

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
#' T <- (numeratorT / exp(1 * X[ ,1, drop = FALSE] + (-0.5 - 1 * X[ ,2, drop = FALSE]) * W)) ^ 2
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
      t0 <- object$t0
    }
    surf1 <- predict(object$fit1, newdata, failure.times = t0)$predictions
    surf0 <- predict(object$fit0, newdata, failure.times = t0)$predictions
    return(surf1 - surf0)
  }
}
