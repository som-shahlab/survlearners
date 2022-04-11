#' @title T-learner of Cox PH
#'
#' @description  T-learner, implemented via Cox proportional hazard models
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
#' @param newX The test data set (covariates only)
#' @examples
#' \donttest{
#' n = 1000; p = 25
#' times = 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[,1] + (-0.5 - 1 * X[,2]) * W))^2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC/(4^2))^(1/2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv_tl_coxph_fit = surv_tl_coxph(X, W, Y, D, times)
#' cate = predict(surv_tl_coxph_fit)
#' cate.test = predict(surv_tl_coxph_fit, X.test)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_tl_coxph <- function(X, W, Y, D, times){

  traindat <- data.frame(Y = Y, D = D, W = W, X)
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]

  # Model for W = 1
  coxph_fit1 <- survival::coxph(survival::Surv(Y, D) ~., data = traindat1)
  bh_dat1 <- survival::basehaz(coxph_fit1, centered = FALSE)
  index <- findInterval(times, bh_dat1$time)
  bh <- bh_dat1[index, 1]
  est_r1 <- predict(coxph_fit1, newdata = data.frame(X), type="risk")
  surf1 <- exp(-bh)^est_r1

  # Model for W = 0
  coxph_fit0 <- survival::coxph(survival::Surv(Y, D) ~., data = traindat0)
  bh_dat0 <- survival::basehaz(coxph_fit0, centered = FALSE)
  index <- findInterval(times, bh_dat0$time)
  bh <- bh_dat0[index, 1]
  est_r0 <- predict(coxph_fit0, newdata = data.frame(X), type="risk")
  surf0 <- exp(-bh)^est_r0

  pred_T_coxph <- surf1 - surf0

  ret <- list(fit1 = coxph_fit1,
              fit0 = coxph_fit0,
              bh1 = bh_dat1,
              bh0 = bh_dat0,
              tau = pred_T_coxph,
              times = times)
  class(ret) <- 'surv_tl_coxph'
  ret
}

#' predict for surv_tl_coxph
#'
#' get estimated tau(X) using the trained surv_tl_coxph model
#'
#' @param object An surv_tl_coxph object
#' @param newx Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param times The prediction time of interest
#' @param ... Additional arguments (currently not used)
#'
#' @examples
#' \donttest{
#' n = 1000; p = 25
#' times = 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[,1] + (-0.5 - 1 * X[,2]) * W))^2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC/(4^2))^(1/2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv_tl_coxph_fit = surv_tl_coxph(X, W, Y, D, times)
#' cate = predict(surv_tl_coxph_fit)
#' cate.test = predict(surv_tl_coxph_fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_tl_coxph = function(object,
                                 newx = NULL,
                                 times = NULL,
                                 ...) {
  if(is.null(newx)){
    return(object$tau)
  }else{
    if(is.null(times)){
      index1 <- findInterval(object$times, object$bh1$time)
      index0 <- findInterval(object$times, object$bh0$time)
    }else{
      index1 <- findInterval(times, object$bh1$time)
      index0 <- findInterval(times, object$bh0$time)
    }

    bh1 <- object$bh1[index1, 1]
    bh0 <- object$bh0[index0, 1]

    est_r1 <- predict(object$fit1, newdata = data.frame(newx), type="risk")
    surf1 <- exp(-bh1)^est_r1

    est_r0 <- predict(object$fit0, newdata = data.frame(newx), type="risk")
    surf0 <- exp(-bh0)^est_r0

    return(surf1 - surf0)
  }
}
