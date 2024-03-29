#' @title S-learner with Lasso
#'
#' @description Estimating conditional average treatment effects (CATEs) for
#' survival outcomes using S-learner with penalized regression models Lasso
#' (implemented via the glmnet package).
#' The CATE is defined as tau(X) = p(Y(1) > t0 | X = x) - p(Y(0) > t0 | X = x),
#' where Y(1) and Y(0) are counterfactual survival times under the treated and controlled arms, respectively.
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param alpha Mix tuning parameter for the elastic net
#' @param k.folds Number of folds for cross validation
#' @param foldid User-supplied foldid. Must have length equal to length(W). If provided, it overrides the k.folds option.
#' @param lambda User-supplied lambda sequence for cross validation
#' @param lambda.choice How to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty.factor User-supplied penalty factor, must be of length the same as number of features in X
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
#' Y <- as.numeric(pmin(failure.time, censor.time))
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.sl.lasso.fit <- surv_sl_lasso(X, Y, W, D, t0)
#' cate <- predict(surv.sl.lasso.fit)
#' cate.test <- predict(surv.sl.lasso.fit, X.test)
#' }
#' @return a surv_sl_lasso object
#' @export
surv_sl_lasso <- function(X, Y, W, D, t0,
                          alpha = 1,
                          k.folds = NULL,
                          foldid = NULL,
                          lambda = NULL,
                          lambda.choice = "lambda.min",
                          penalty.factor = NULL) {

  input <- sanitize_input(X, Y, W, D)
  X <- input$X
  W <- input$W
  Y <- input$Y
  D <- input$D

  x.mean <- colMeans(X)
  x.sd <- apply(X, 2, sd)
  x.sd[x.sd == 0] <- 1 # in case Xj is constant.
  x.scl <- scale(X, center = x.mean, scale = x.sd)

  nobs <- nrow(x.scl)
  pobs <- ncol(x.scl)

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

  x.scl.tilde <- cbind(as.numeric(W - 0.5) * cbind(1, x.scl), x.scl)
  x.scl.pred1 <- cbind(0.5 * cbind(1, x.scl), x.scl)
  x.scl.pred0 <- cbind(-0.5 * cbind(1, x.scl), x.scl)

  if (is.null(penalty.factor) || (length(penalty.factor) != pobs)) {
    if (!is.null(penalty.factor) && length(penalty.factor) != 2 * pobs + 1) {
      warning("penalty.factor supplied is not 1 plus 2 times the number of columns in X. Using all ones instead.")
    }
    penalty.factor <- c(0, rep(1, 2 * pobs))
  }
  x.scl.tilde <- as.matrix(data.frame(x.scl.tilde))
  tau.fit <- glmnet::cv.glmnet(x.scl.tilde,
                               survival::Surv(Y, D),
                               family = "cox",
                               foldid = foldid,
                               lambda = lambda,
                               penalty.factor = penalty.factor,
                               standardize = FALSE,
                               alpha = alpha)

  s.beta <- as.vector(coef(tau.fit, s = lambda.choice))
  link1 <- exp(x.scl.pred1 %*% s.beta)
  link0 <- exp(x.scl.pred0 %*% s.beta)
  S0.t <- base_surv(fit = tau.fit,
                    Y = Y,
                    D = D,
                    X = x.scl.tilde,
                    lambda = tau.fit$lambda.min)
  index <- findInterval(t0, S0.t$time)
  if (index == 0) {
    S0 <- 1
  } else {
    S0 <- S0.t[index, ]$survival
  }

  surv1 <- S0 ^ link1
  surv0 <- S0 ^ link0

  tau.hat <- as.numeric(surv1 - surv0)

  ret <- list(tau.fit = tau.fit,
              x.org = x.scl.tilde,
              y.org = Y,
              D.org = D,
              s.beta = s.beta,
              S0.t = S0.t,
              t0 = t0,
              tau.hat = tau.hat,
              lambda.choice = lambda.choice,
              standardization = list(x.mean = x.mean, x.sd = x.sd))

  class(ret) <- "surv_sl_lasso"
  ret
}

#' Predict with a S-learner with Lasso
#'
#' Obtain estimated tau(X) using a trained S-learner with Lasso model
#'
#' Remark: CATE predictions can be made at any time point on the estimated survival curve
#'
#' @param object A surv_sl_lasso object
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
#' Y <- as.numeric(pmin(failure.time, censor.time))
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.sl.lasso.fit <- surv_sl_lasso(X, Y, W, D, t0)
#' cate <- predict(surv.sl.lasso.fit)
#' cate.test <- predict(surv.sl.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_sl_lasso <- function(object,
                                  newdata = NULL,
                                  t0 = NULL,
                                  ...) {
  if (!is.null(newdata)) {
    newdata <- sanitize_x(newdata)
    newdata.scl <- scale(newdata, center = object$standardization$x.mean, scale = object$standardization$x.sd)
    newdata.scl.pred1 <- cbind(0.5 * cbind(1, newdata.scl), newdata.scl)
    newdata.scl.pred0 <- cbind(-0.5 * cbind(1, newdata.scl), newdata.scl)

    link1 <- exp(newdata.scl.pred1 %*% object$s.beta)
    link0 <- exp(newdata.scl.pred0 %*% object$s.beta)

    if (is.null(t0)) {
    t0 <- object$t0
    }
    index <- findInterval(t0, object$S0.t$time)
    if (index == 0) {
      S0 <- 1
    } else {
      S0 <- object$S0.t[index, ]$survival
    }

    surv1 <- S0^link1
    surv0 <- S0^link0

    tau.hat <- as.numeric(surv1 - surv0)
  }
  else {
    tau.hat <- object$tau.hat
  }
  return(tau.hat)
}
