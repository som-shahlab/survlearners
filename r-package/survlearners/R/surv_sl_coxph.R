#' @title S-learner with Cox regression
#'
#' @description Estimating conditional average treatment effects (CATEs) for
#' survival outcomes using S-learner with Cox proportional hazard models
#' (implemented via the survival package).
#' The CATE is defined as tau(X) = p(Y(1) > t0 | X = x) - p(Y(0) > t0 | X = x),
#' where Y(1) and Y(0) are counterfactual survival times under the treated and controlled arms, respectively.
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param weights The (optional) case weights
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
#' surv.sl.coxph.fit <- surv_sl_coxph(X, Y, W, D, t0)
#' cate <- predict(surv.sl.coxph.fit)
#' cate.test <- predict(surv.sl.coxph.fit, X.test)
#' }
#' @return A surv_sl_coxph object
#' @export
surv_sl_coxph <- function(X, Y, W, D, t0, weights = NULL) {

  input <- sanitize_input(X, Y, W, D)
  X <- input$X
  W <- input$W
  Y <- input$Y
  D <- input$D

  x.tilde <- data.frame(as.numeric(W - 0.5) * cbind(1, X), X)
  x.pred1 <- data.frame(0.5 * cbind(1, X), X)
  x.pred0 <- data.frame(-0.5 * cbind(1, X), X)

  colnames(x.tilde) <- colnames(x.pred1) <- colnames(x.pred0) <- paste0("v", 1:dim(x.tilde)[2])
  formula <- as.formula(paste0("survival::Surv(Y, D) ~ ", paste(colnames(x.tilde), sep=" ", collapse = "+")))
  tmpdat <- data.frame(Y, D, x.tilde)

  if(is.null(weights)) {
    tau.fit <- survival::coxph(formula, data = tmpdat)
  } else {
    # add simple check length(weights) == length(Y)
    tau.fit <- survival::coxph(formula, weights = weights, data = tmpdat)
  }

  bh.dat <- survival::basehaz(tau.fit, centered = FALSE)
  index <- findInterval(t0, bh.dat$time)
  if (index == 0) {
    bh <- 0
  } else {
    bh <- bh.dat[index, 1]
  }

  link1 <- exp(as.matrix(x.pred1) %*% tau.fit$coefficients)
  link0 <- exp(as.matrix(x.pred0) %*% tau.fit$coefficients)

  est.S1.cvd <- exp(-bh) ^ link1
  est.S0.cvd <- exp(-bh) ^ link0

  tau.hat <- est.S1.cvd - est.S0.cvd

  ret <- list(tau.fit = tau.fit,
              bh = bh,
              s.beta = tau.fit$coefficients,
              t0 = t0,
              tau.hat = tau.hat)

  class(ret) <- "surv_sl_coxph"
  ret
}

#' Predict for a S-learner with Cox regression
#'
#' Obtain estimated tau(X) using a trained S-learner with Cox regression
#'
#' Remark: CATE predictions can be made at any time point on the estimated survival curve
#'
#' @param object A surv_sl_coxph object
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
#' surv.sl.coxph.fit <- surv_sl_coxph(X, Y, W, D, t0)
#' cate <- predict(surv.sl.coxph.fit)
#' cate.test <- predict(surv.sl.coxph.fit, X.test)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_sl_coxph <- function(object,
                                  newdata = NULL,
                                  t0 = NULL,
                                  ...) {
  if (!is.null(newdata)) {
    newdata <- sanitize_x(newdata)

    if (is.null(t0)) {
    t0 <- object$t0
    }

    bh.dat <- survival::basehaz(object$tau.fit, centered = FALSE)
    index <- findInterval(t0, bh.dat$time)
    if (index == 0) {
      bh <- 0
    } else {
      bh <- bh.dat[index, 1]
    }

    x.pred1 <- data.frame(0.5 * cbind(1, newdata), newdata)
    x.pred0 <- data.frame(-0.5 * cbind(1, newdata), newdata)

    link1 <- exp(as.matrix(x.pred1) %*% object$s.beta)
    link0 <- exp(as.matrix(x.pred0) %*% object$s.beta)

    est.S1.cvd <- exp(-bh) ^ link1
    est.S0.cvd <- exp(-bh) ^ link0

    tau.hat <- est.S1.cvd - est.S0.cvd
  }
  else {
    tau.hat <- object$tau.hat
  }
  return(tau.hat)
}
