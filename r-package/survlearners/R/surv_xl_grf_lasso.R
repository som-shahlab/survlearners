#' @title X-learner with random survival forest and Lasso
#'
#' @description Estimating conditional average treatment effects (CATEs) for
#' survival outcomes using X-learner with random survival for nuisance outcomes
#' (implemented via the grf package) and penalized regression Lasso for CATEs
#' (implemented via the glmnet package).
#' The CATE is defined as tau(X) = p(Y(1) > t0 | X = x) - p(Y(0) > t0 | X = x),
#' where Y(1) and Y(0) are counterfactual survival times under the treated and controlled arms, respectively.
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param W.hat The propensity score
#' @param cen.fit The choice of model fitting for censoring
#' @param k.folds The number of folds for estimating nuisance parameters via cross-fitting
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
#' surv.xl.grf.lasso.fit <- surv_xl_grf_lasso(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.xl.grf.lasso.fit)
#' cate.test <- predict(surv.xl.grf.lasso.fit, X.test)
#' }
#' @return A surv_xl_grf_lasso object
#' @export
surv_xl_grf_lasso <- function(X, Y, W, D, t0, W.hat = NULL, cen.fit = "Kaplan-Meier",
                              k.folds = 10, new.args.grf.nuisance = list()) {

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

  # fit model on W == 1
  grffit1 <- do.call(grf::survival_forest, c(list(X = X[W == 1,, drop = FALSE], Y = Y[W == 1], D = D[W == 1]), args.grf.nuisance))
  surf1 <- predict(grffit1, X, failure.times = t0)$predictions

  # fit model on W == 0
  grffit0 <- do.call(grf::survival_forest, c(list(X = X[W == 0,, drop = FALSE], Y = Y[W == 0], D = D[W == 0]), args.grf.nuisance))
  surf0 <- predict(grffit0, X, failure.times = t0)$predictions

  Tgrf1 <- 1 - surf1
  Tgrf0 <- 1 - surf0

  # IPCW weights
  U <- pmin(Y, t0)                         # truncated follow-up time by t0
  if (cen.fit == "Kaplan-Meier") {
    fold.id <- sample(rep(seq(k.folds), length = nrow(X)))
    C.hat <- rep(NA, length(fold.id))
    for (z in 1:k.folds) {
      c.fit <- survival::survfit(survival::Surv(Y[!fold.id == z], 1 - D[!fold.id == z]) ~ 1)
      kmc <- summary(c.fit, times = U[fold.id == z])
      C.hat[fold.id == z] <- kmc$surv[match(U[fold.id == z], kmc$time)]
    }
  } else if (cen.fit == "survival.forest") {
    args.grf.nuisance$alpha <- 0.05
    c.fit <- do.call(grf::survival_forest, c(list(X = cbind(W, X), Y = Y, D = 1 - D), args.grf.nuisance))
    C.hat <- predict(c.fit, failure.times = U, prediction.times = "time")$predictions
  }
  if (any(C.hat == 0)) {
    stop("Some or all uncensored probabilities are exactly zeros. Check input variables or consider adjust the time of interest t0.")
  }

  # Propensity score
  if (is.null(W.hat)) {
    stop("propensity score needs to be supplied")
  } else if (length(W.hat) == 1) {
    W.hat <- rep(W.hat, length(W))
  } else if (length(W.hat) != nrow(X)) {
    stop("W.hat has incorrect length.")
  }

  # CATE function
   D.t0 <- D
   D.t0[D == 1 & Y > t0] <- 0
   D.t0 <- D.t0[D == 1 | Y > t0]
   W.t0 <- W[D == 1 | Y > t0]
   X.t0 <- X[D == 1 | Y > t0,, drop = FALSE]
   Tgrf0.t0 <- Tgrf0[D == 1 | Y > t0]
   Tgrf1.t0 <- Tgrf1[D == 1 | Y > t0]
   sample.weights.t0 <- 1 / (C.hat * W.hat)[D == 1 | Y > t0]

  foldid <- sample(rep(seq(k.folds), length = length(D.t0[W.t0 == 1])))
  tau.fit1 <- glmnet::cv.glmnet(X.t0[W.t0 == 1,, drop = FALSE],
                                D.t0[W.t0 == 1] - Tgrf0.t0[W.t0 == 1],
                                weights = sample.weights.t0[W.t0 == 1],
                                foldid = foldid,
                                alpha = 1)
  XLtau1 <- as.vector(-predict(tau.fit1, X, s = "lambda.min"))

  foldid <- sample(rep(seq(k.folds), length = length(D.t0[W.t0 == 0])))
  tau.fit0 <- glmnet::cv.glmnet(X.t0[W.t0 == 0,, drop = FALSE],
                                Tgrf1.t0[W.t0 == 0] - D.t0[W.t0 == 0],
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
  class(ret) <- "surv_xl_grf_lasso"
  ret
}

#' Predict with a X-learner with random survival forest and Lasso
#'
#' Obtain estimated tau(X) using a trained X-learner with random survival forest and Lasso model
#'
#' Remark: CATE predictions can only be made at the time point used to define the outcome in the trained model
#'
#' @param object An surv_xl_grf_lasso object
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
#' T <- (numeratorT / exp(1 * X[ ,1, drop = FALSE] + (-0.5 - 1 * X[ ,2, drop = FALSE]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.xl.grf.lasso.fit <- surv_xl_grf_lasso(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.xl.grf.lasso.fit)
#' cate.test <- predict(surv.xl.grf.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_xl_grf_lasso <- function(object,
                                      newdata = NULL,
                                      ...) {
  if (is.null(newdata)) {
    return(object$tau.hat)
  } else {
    XLtau1 <- as.vector(-predict(object$tau.fit1, newdata, s = "lambda.min"))
    XLtau0 <- as.vector(-predict(object$tau.fit0, newdata, s = "lambda.min"))
    W.hat <- rep(object$W.hat[1], nrow(newdata))
    return(as.vector(XLtau1 * (1 - W.hat) + XLtau0 * W.hat))
  }
}
