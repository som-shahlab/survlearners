#' @title F-learner of grf
#'
#' @description  F-learner, implemented via survival_forest in the grf package
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param W.hat The propensity score
#' @param cen.fit The choice of model fitting for censoring
#' @param k.folds The number of folds for estimating nuisance parameters via cross-fitting
#' @param args.grf.nuisance Input arguments for a grf model that estimates nuisance parameters
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
#' surv.fl.grf.fit <- surv_fl_grf(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.fl.grf.fit)
#' cate.test <- predict(surv.fl.grf.fit, X.test)
#' }
#' @return A surv_fl_grf object
#' @export
surv_fl_grf <- function(X, Y, W, D, t0, W.hat = NULL, cen.fit = "Kaplan-Meier",
                        k.folds = 10, args.grf.nuisance = list()) {

  # IPCW weights
  Q <- as.numeric(D == 1 | Y > t0)         # indicator for uncensored at t0
  U <- pmin(Y, t0)                         # truncated follow-up time by t0
  if (cen.fit == "Kaplan-Meier") {
    fold.id <- sample(rep(seq(k.folds), length = nrow(X)))
    C.hat <- rep(NA, length(fold.id))
    for (z in 1:k.folds) {
      c.fit <- survival::survfit(survival::Surv(Y[!fold.id == z], 1 - Q[!fold.id == z]) ~ 1)
      C.hat[fold.id == z] <- summary(c.fit, times = U[fold.id == z])$surv
    }
  } else if (cen.fit == "survival.forest") {
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
    c.fit <- do.call(grf::survival_forest, c(list(X = cbind(W, X), Y = Y, D = 1 - Q), args.grf.nuisance))
    C.hat <- predict(c.fit)$predictions
    cen.times.index <- findInterval(U, c.fit$failure.times)
    C.hat <- C.hat[cbind(1:length(U), cen.times.index)]
  }

  # Propensity score
  if (is.null(W.hat)) {
    stop("propensity score needs to be supplied")
  } else {
    W.hat <- rep(W.hat, length(Y))
  }

  # CATE function
  D.t0 <- D
  D.t0[D == 1 & Y > t0] <- 0
  D.t0 <- D.t0[D == 1 | Y > t0]
  W.t0 <- W[D == 1 | Y > t0]
  X.t0 <- X[D == 1 | Y > t0,, drop = FALSE]
  sample.weights.t0 <- 1 / C.hat[D == 1 | Y > t0]
  W.hat.t0 <- W.hat[D == 1 | Y > t0]

  Z <- W.t0 * D.t0 / W.hat.t0 - (1 - W.t0) * D.t0 / (1 - W.hat.t0)
  tau.fit <- grf::regression_forest(X.t0, Z, sample.weights = sample.weights.t0)
  tau.hat <- -predict(tau.fit, X.t0)

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
#' surv.fl.grf.fit <- surv_fl_grf(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.fl.grf.fit)
#' cate.test <- predict(surv.fl.grf.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_fl_grf <- function(object,
                                newdata = NULL,
                                ...) {
  if (is.null(newdata)) {
    return(object$tau.hat)
  } else {
    return(-predict(object$tau.fit, newdata))
  }
}
