#' @title T-learner with Lasso
#'
#' @description Estimating conditional average treatment effects (CATEs) for
#' survival outcomes using T-learner with penalized regression models Lasso
#' (implemented via the glmnet package).
#' The CATE is defined as tau(X) = p(Y(1) > t0 | X = x) - p(Y(0) > t0 | X = x),
#' where Y(1) and Y(0) are counterfactual survival times under the treated and controlled arms, respectively.
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param k.folds The number of folds for estimating nuisance parameters via cross-fitting
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
#' surv.tl.lasso.fit <- surv_tl_lasso(X, Y, W, D, t0)
#' cate <- predict(surv.tl.lasso.fit)
#' cate.test <- predict(surv.tl.lasso.fit, X.test)
#' }
#' @return A surv_tl_lasso object
#' @export
surv_tl_lasso <- function(X, Y, W, D, t0, k.folds = 10) {

  # Model for W = 1
  foldid <- sample(rep(seq(k.folds), length = length(Y[W == 1])))
  x1 <- as.matrix(data.frame(X[W == 1,, drop = FALSE]))
  lasso.fit1 <- glmnet::cv.glmnet(x1,
                                  survival::Surv(Y[W == 1], D[W == 1]),
                                  family = "cox",
                                  alpha = 1,
                                  foldid = foldid)

  bsurv1 <- base_surv(fit = lasso.fit1,
                      Y = Y[W == 1],
                      D = D[W == 1],
                      X = x1,
                      lambda = lasso.fit1$lambda.min)

  surf1 <- pred_surv(fit = lasso.fit1,
                     S0 = bsurv1,
                     X = X,
                     t0 = t0,
                     lambda = lasso.fit1$lambda.min)


  # Model for W = 0
  foldid <- sample(rep(seq(k.folds), length = length(Y[W == 0])))
  x0 <- as.matrix(data.frame(X[W == 0,, drop = FALSE]))
  lasso.fit0 <- glmnet::cv.glmnet(x0,
                                  survival::Surv(Y[W == 0], D[W == 0]),
                                  family = "cox",
                                  alpha = 1,
                                  foldid = foldid)

  bsurv0 <- base_surv(fit = lasso.fit0,
                      Y = Y[W == 0],
                      D = D[W == 0],
                      X = x0,
                      lambda = lasso.fit0$lambda.min)

  surf0 <- pred_surv(fit = lasso.fit0,
                     S0 = bsurv0,
                     X = X,
                     t0 = t0,
                     lambda = lasso.fit0$lambda.min)

  tau.hat <- surf1 - surf0

  ret <- list(fit1 = lasso.fit1,
              fit0 = lasso.fit0,
              bsurv1 = bsurv1,
              bsurv0 = bsurv0,
              tau.hat = tau.hat,
              t0 = t0)
  class(ret) <- "surv_tl_lasso"
  ret
}

#' Predict with a T-learner with Lasso
#'
#' Obtain estimated tau(X) using a trained T-learner with Lasso model
#'
#' Remark: CATE predictions can be made at any time point on the estimated survival curve
#'
#' @param object An surv_tl_lasso object
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
#' surv.tl.lasso.fit <- surv_tl_lasso(X, Y, W, D, t0)
#' cate <- predict(surv.tl.lasso.fit)
#' cate.test <- predict(surv.tl.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_tl_lasso <- function(object,
                                  newdata = NULL,
                                  t0 = NULL,
                                  ...) {
  if (is.null(newdata)) {
    return(object$tau.hat)
  } else {
    if (is.null(t0)) {
    t0 <- object$t0
    }
    surf1 <- pred_surv(fit = object$fit1,
                       S0 = object$bsurv1,
                       X = newdata,
                       t0 = t0,
                       lambda = object$fit1$lambda.min)

    surf0 <- pred_surv(fit = object$fit0,
                       S0 = object$bsurv0,
                       X = newdata,
                       t0 = t0,
                       lambda = object$fit0$lambda.min)
    return(surf1 - surf0)
  }
}
