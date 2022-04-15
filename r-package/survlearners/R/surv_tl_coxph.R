#' @title T-learner of Cox PH
#'
#' @description  T-learner, implemented via Cox proportional hazard models
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
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
#' surv.tl.coxph.fit <- surv_tl_coxph(X, Y, W, D, t0)
#' cate <- predict(surv.tl.coxph.fit)
#' cate.test <- predict(surv.tl.coxph.fit, X.test)
#' }
#' @return A surv_tl_coxph object
#' @export
surv_tl_coxph <- function(X, Y, W, D, t0) {

  traindat <- data.frame(Y = Y, D = D, W = W, X)
  traindat1 <- traindat[traindat$W == 1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W == 0, !colnames(traindat) %in% c("W")]

  # Model for W = 1
  coxph.fit1 <- survival::coxph(survival::Surv(Y, D) ~ ., data = traindat1)
  bh.dat1 <- survival::basehaz(coxph.fit1, centered = FALSE)
  index <- findInterval(t0, bh.dat1$time)
  bh <- bh.dat1[index, 1]
  est.r1 <- predict(coxph.fit1, newdata = data.frame(X), type="risk")
  surf1 <- exp(-bh) ^ est.r1

  # Model for W = 0
  coxph.fit0 <- survival::coxph(survival::Surv(Y, D) ~ ., data = traindat0)
  bh.dat0 <- survival::basehaz(coxph.fit0, centered = FALSE)
  index <- findInterval(t0, bh.dat0$time)
  bh <- bh.dat0[index, 1]
  est.r0 <- predict(coxph.fit0, newdata = data.frame(X), type="risk")
  surf0 <- exp(-bh) ^ est.r0

  tau.hat <- surf1 - surf0

  ret <- list(fit1 = coxph.fit1,
              fit0 = coxph.fit0,
              bh1 = bh.dat1,
              bh0 = bh.dat0,
              tau.hat = tau.hat,
              t0 = t0)
  class(ret) <- "surv_tl_coxph"
  ret
}

#' predict for surv_tl_coxph
#'
#' get estimated tau(X) using the trained surv_tl_coxph model
#'
#' @param object An surv_tl_coxph object
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
#' surv.tl.coxph.fit <- surv_tl_coxph(X, Y, W, D, t0)
#' cate <- predict(surv.tl.coxph.fit)
#' cate.test <- predict(surv.tl.coxph.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_tl_coxph <- function(object,
                                  newdata = NULL,
                                  t0 = NULL,
                                  ...) {
  if (is.null(newdata)) {
    return(object$tau.hat)
  } else {
    if (is.null(t0)) {
      index1 <- findInterval(object$t0, object$bh1$time)
      index0 <- findInterval(object$t0, object$bh0$time)
    } else {
      index1 <- findInterval(t0, object$bh1$time)
      index0 <- findInterval(t0, object$bh0$time)
    }

    bh1 <- object$bh1[index1, 1]
    bh0 <- object$bh0[index0, 1]

    est.r1 <- predict(object$fit1, newdata = data.frame(newdata), type="risk")
    surf1 <- exp(-bh1) ^ est.r1

    est.r0 <- predict(object$fit0, newdata = data.frame(newdata), type="risk")
    surf0 <- exp(-bh0) ^ est.r0

    return(surf1 - surf0)
  }
}
