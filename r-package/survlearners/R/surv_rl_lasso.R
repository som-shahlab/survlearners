#' @title R-learner with Lasso
#'
#' @description Estimating conditional average treatment effects (CATEs) for
#' survival outcomes using R-learner with penalized regression models Lasso
#' (implemented via the glmnet package).
#' The CATE is defined as tau(X) = p(Y(1) > t0 | X = x) - p(Y(0) > t0 | X = x),
#' where Y(1) and Y(0) are counterfactual survival times under the treated and controlled arms, respectively.
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param k.folds Number of folds for cross validation
#' @param foldid User-supplied foldid. Must have length equal to length(W). If provided, it overrides the k.folds option.
#' @param W.hat Propensity score
#' @param Y.hat Conditional mean outcome E(Y|X)
#' @param C.hat Censoring weights
#' @param lambda.choice How to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param new.args.lasso.nuisance Input arguments for a lasso model that estimates nuisance parameters
#' @param new.args.grf.nuisance Input arguments for a grf model that estimates nuisance parameters
#' @param new.args.lasso.tau Input arguments for a lasso model that estimates CATE
#' @param cen.fit The choice of model fitting for censoring
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
#' surv.rl.lasso.fit <- surv_rl_lasso(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.rl.lasso.fit)
#' cate.test <- predict(surv.rl.lasso.fit, X.test)
#' }
#' @return A surv_rl_lasso object
#' @export
surv_rl_lasso <- function(X, Y, W, D,
                          t0 = NULL,
                          k.folds = 10,
                          foldid = NULL,
                          W.hat = NULL,
                          Y.hat = NULL,
                          C.hat = NULL,
                          lambda.choice = "lambda.min",
                          new.args.lasso.nuisance = list(),
                          new.args.grf.nuisance = list(),
                          new.args.lasso.tau = list(),
                          cen.fit = "Kaplan-Meier") {

    input <- sanitize_input(X, Y, W, D)
    X <- input$X
    W <- input$W
    Y <- input$Y
    D <- input$D

    if (is.null(foldid) || length(foldid) != length(W)) {

      if (!is.null(foldid) && length(foldid) != length(W)) {
        warning("supplied foldid does not have the same length ")
      }

      if (is.null(k.folds)) {
          k.folds <- floor(max(3, min(10, length(W) / 4)))
      }

      # fold ID for cross-validation; balance treatment assignments
      foldid <- sample(rep(seq(k.folds), length = length(W)))

    }

    if (is.null(W.hat)) {
      stop("propensity score needs to be supplied")
    } else if (length(W.hat) == 1) {
      W.hat <- rep(W.hat, nrow(X))
    } else if (length(W.hat) != nrow(X)) {
      stop("W.hat has incorrect length.")
    }

    args.lasso.nuisance <- list(family = "cox",
                                nfolds = k.folds,
                                alpha = 1,
                                lambda = NULL)
    args.lasso.nuisance[names(new.args.lasso.nuisance)] <- new.args.lasso.nuisance

    if (is.null(Y.hat)) {
    foldid <- sample(rep(seq(k.folds), length = length(W)))
    survt1 <- survt0 <- rep(NA, length(W))
    for (k in 1:k.folds) {
      X1 <- as.matrix(data.frame(X[(!foldid == k) & W == 1,, drop = FALSE]))
      y1 <- Y[(!foldid == k) & W == 1]
      D1 <- D[(!foldid == k) & W == 1]
      y1.fit <- do.call(glmnet::cv.glmnet, c(list(x = X1, y = survival::Surv(y1, D1)), args.lasso.nuisance))
      S0.1 <- base_surv(y1.fit, y1, D1, X1, lambda = y1.fit$lambda.min)
      survt1[foldid == k] <- pred_surv(y1.fit, S0.1, X[foldid == k, ], t0 = t0, lambda = y1.fit$lambda.min)

      X0 <- as.matrix(data.frame(X[(!foldid == k) & W == 0,, drop = FALSE]))
      y0 <- Y[(!foldid == k) & W == 0]
      D0 <- D[(!foldid == k) & W == 0]
      y0.fit <- do.call(glmnet::cv.glmnet, c(list(x = X0, y = survival::Surv(y0, D0)), args.lasso.nuisance))
      S0.0 <- base_surv(y0.fit, y0, D0, X0, lambda = y0.fit$lambda.min)
      survt0[foldid == k] <- pred_surv(y0.fit, S0.0, X[foldid == k, ], t0 = t0, lambda = y0.fit$lambda.min)
    }
    Y.hat  <- W.hat * survt1 + (1 - W.hat) * survt0
    } else {
      y1.fit <- NULL
      y0.fit <- NULL
    }

    args.grf.nuisance <- list(failure.times = NULL,
                              num.trees = max(50, 2000 / 4),
                              min.node.size = 15,
                              honesty = TRUE,
                              honesty.fraction = 0.5,
                              honesty.prune.leaves = TRUE,
                              alpha = 0.05,
                              prediction.type = "Nelson-Aalen",
                              compute.oob.predictions = TRUE)
    args.grf.nuisance[names(new.args.grf.nuisance)] <- new.args.grf.nuisance

    if (is.null(C.hat)) {
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
        c.fit <- do.call(grf::survival_forest, c(list(X = cbind(W, X), Y = Y, D = 1 - D), args.grf.nuisance))
        C.hat <- predict(c.fit, failure.times = U, prediction.times = "time")$predictions
      }
    } else {
      c.fit <- NULL
    }
    if (any(C.hat == 0)) {
      stop("Some or all uncensored probabilities are exactly zeros. Check input variables or consider adjust the time of interest t0.")
    }

    # CATE function
    D.t0 <- D
    D.t0[D == 1 & Y > t0] <- 0
    D.t0 <- D.t0[D == 1 | Y > t0]
    Y.hat.t0 <- Y.hat[D == 1 | Y > t0]
    W.t0 <- W[D == 1 | Y > t0]
    W.hat.t0 <- W.hat[D == 1 | Y > t0]
    X.t0 <- X[D == 1 | Y > t0,, drop = FALSE]
    sample.weights.t0 <- (W.t0 - W.hat.t0)^2 / C.hat[D == 1 | Y > t0]
    y.tilde <- ((1 - D.t0) - Y.hat.t0) / (W.t0 - W.hat.t0)
    foldid <- sample(rep(seq(k.folds), length = nrow(X.t0)))
    args.lasso.tau <- list(weights = sample.weights.t0,
                           foldid = foldid,
                           alpha = 1,
                           lambda = NULL)
    args.lasso.tau[names(new.args.lasso.nuisance)] <- new.args.lasso.tau
    tau.fit <- do.call(glmnet::cv.glmnet, c(list(x = X.t0, y = y.tilde), args.lasso.tau))
    tau.hat <- predict(tau.fit, newx = X, s = lambda.choice)

    ret <- list(tau.fit = tau.fit,
                y1.fit = y1.fit,
                y0.fit = y0.fit,
                c.fit = c.fit,
                W.hat = W.hat,
                Y.hat = Y.hat,
                C.hat = C.hat,
                tau.hat = tau.hat)
    class(ret) <- "surv_rl_lasso"
    ret
}


#' Predict with a R-learner with Lasso
#'
#' Obtain estimated tau(X) using a trained R-learner with Lasso model
#'
#' Remark: CATE predictions can only be made at the time point used to define the outcome in the trained model
#'
#' @param object An surv_rl_lasso object
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
#' surv.rl.lasso.fit <- surv_rl_lasso(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.rl.lasso.fit)
#' cate.test <- predict(surv.rl.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_rl_lasso <- function(object,
                                  newdata = NULL,
                                  ...) {
  if (!is.null(newdata)) {
    newdata <- sanitize_x(newdata)
    tau.hat <- predict(object$tau.fit, newx = newdata, s = "lambda.min")
  } else {
    tau.hat <- object$tau.hat
  }
  return(tau.hat)
}
