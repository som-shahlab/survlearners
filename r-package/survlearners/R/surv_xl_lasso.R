#' @title X-learner of lasso
#'
#' @description  X-learner, implemented via glmnet (lasso) with "coxph" distribution
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
#' T <- (numeratorT / exp(1 * X[ ,1, drop = FALSE] + (-0.5 - 1 * X[ ,2, drop = FALSE]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.xl.lasso.fit <- surv_xl_lasso(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.xl.lasso.fit)
#' cate.test <- predict(surv.xl.lasso.fit, X.test)
#' }
#' @return A surv_xl_lasso object
#' @export
surv_xl_lasso <- function(X, Y, W, D, t0, W.hat = NULL, cen.fit = "Kaplan-Meier",
                          k.folds = 10, args.grf.nuisance = list()) {

  # fit model on W == 1 (cross-fitted using "preval" in glmnet)
  foldid1 <- sample(rep(seq(k.folds), length = length(Y[W == 1])))
  x1 <- as.matrix(data.frame(X[W == 1,, drop = FALSE]))
  lasso.fit1 <- glmnet::cv.glmnet(x1,
                                  survival::Surv(Y[W == 1], D[W == 1]),
                                  family = "cox",
                                  alpha = 1,
                                  foldid = foldid1)
  S0 <- base_surv(lasso.fit1, Y[W == 1], D[W == 1], x1, lambda = lasso.fit1$lambda.min)
  surf1 <- pred_surv(fit = lasso.fit1, S0 = S0, X = X, t0 = t0, lambda = lasso.fit1$lambda.min)

  # fit model on W == 0 (cross-fitted using "preval" in glmnet)
  foldid0 <- sample(rep(seq(k.folds), length = length(Y[W == 0])))
  x0 <- as.matrix(data.frame(X[W == 0,, drop = FALSE]))
  lasso.fit0 <- glmnet::cv.glmnet(x0,
                                  survival::Surv(Y[W == 0], D[W == 0]),
                                  family = "cox",
                                  alpha = 1,
                                  foldid = foldid0)
  S0 <- base_surv(lasso.fit0, Y[W == 0], D[W == 0], x0, lambda = lasso.fit0$lambda.min)
  surf0 <- pred_surv(fit = lasso.fit0, S0 = S0, X = X, t0 = t0, lambda = lasso.fit0$lambda.min)

  Tlasso1 <- 1 - surf1
  Tlasso0 <- 1 - surf0

  # IPCW weights
  U <- pmin(Y, t0)                         # truncated follow-up time by t0
  if (cen.fit == "Kaplan-Meier") {
    fold.id <- sample(rep(seq(k.folds), length = nrow(X)))
    C.hat <- rep(NA, length(fold.id))
    for (z in 1:k.folds) {
      c.fit <- survival::survfit(survival::Surv(Y[!fold.id == z], 1 - D[!fold.id == z]) ~ 1)
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
    c.fit <- do.call(grf::survival_forest, c(list(X = cbind(W, X), Y = Y, D = 1 - D), args.grf.nuisance))
    C.hat <- predict(c.fit, failure.times = U, prediction.times = "time")$predictions
  }
  if (any(C.hat == 0)) {
    stop("Some or all uncensored probabilities are exactly zeros. Check input variables or consider adjust the time of interest t0.")
  }

  # Propensity score (cross-fitted using "preval" in glmnet)
  if (is.null(W.hat)) {
    stop("propensity score needs to be supplied")
  } else {
    W.hat <- rep(W.hat, length(W))
  }

  # X-learner
   D.t0 <- D
   D.t0[D == 1 & Y > t0] <- 0
   D.t0 <- D.t0[D == 1 | Y > t0]
   W.t0 <- W[D == 1 | Y > t0]
   X.t0 <- X[D == 1 | Y > t0,, drop = FALSE]
   Tlasso0.t0 <- Tlasso0[D == 1 | Y > t0]
   Tlasso1.t0 <- Tlasso1[D == 1 | Y > t0]
   sample.weights.t0 <- 1 / (C.hat * W.hat)[D == 1 | Y > t0]

   foldid <- sample(rep(seq(k.folds), length = length(D.t0[W.t0 == 1])))
   tau.fit1 <- glmnet::cv.glmnet(X.t0[W.t0 == 1,, drop = FALSE],
                                 D.t0[W.t0 == 1] - Tlasso0.t0[W.t0 == 1],
                                 weights = sample.weights.t0[W.t0 == 1],
                                 foldid = foldid,
                                 alpha = 1)
   XLtau1 <- as.vector(-predict(tau.fit1, X, s = "lambda.min"))

   foldid <- sample(rep(seq(k.folds), length = length(D.t0[W.t0 == 0])))
   tau.fit0 <- glmnet::cv.glmnet(X.t0[W.t0 == 0,, drop = FALSE],
                                 Tlasso1.t0[W.t0 == 0] - D.t0[W.t0 == 0],
                                 weights = sample.weights.t0[W.t0 == 0],
                                 foldid = foldid,
                                 alpha = 1)
   XLtau0 <- as.vector(-predict(tau.fit0, X, s = "lambda.min"))

  # weighted CATE
  tau.hat <- as.vector(XLtau1 * (1 - W.hat) + XLtau0 * W.hat)

  ret <- list(tau.fit1 = tau.fit1,
              tau.fit0 = tau.fit0,
              W.hat = W.hat,
              tau.hat = tau.hat)
  class(ret) <- "surv_xl_lasso"
  ret
}

#' predict for surv_xl_lasso
#'
#' get estimated tau(X) using the trained surv_xl_lasso model
#'
#' @param object An surv_xl_lasso object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param W.hat The propensity score
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
#' surv.xl.lasso.fit <- surv_xl_lasso(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.xl.lasso.fit)
#' cate.test <- predict(surv.xl.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_xl_lasso <- function(object,
                                  newdata = NULL,
                                  W.hat = NULL,
                                  ...) {
  if (is.null(newdata)) {
    return(object$tau.hat)
  } else {
    XLtau1 <- as.vector(-predict(object$tau.fit1, newdata, s = "lambda.min"))
    XLtau0 <- as.vector(-predict(object$tau.fit0, newdata, s = "lambda.min"))
    if (is.null(W.hat)) {
      W.hat <- rep(object$W.hat[1], nrow(newdata))
    } else {
      if (length(W.hat) == 1) {
      W.hat <- rep(W.hat, nrow(newdata))
      }
    }
    return(as.vector(XLtau1 * (1 - W.hat) + XLtau0 * W.hat))
  }
}
