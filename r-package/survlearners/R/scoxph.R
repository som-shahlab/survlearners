#' @title S-learner of Cox PH
#'
#' @description  S-learner, implemented via Cox proportional hazard models
#'
#' @param x The baseline covariates
#' @param w The treatment variable (0 or 1)
#' @param y The follow-up time
#' @param D The event indicator
#' @param times The prediction time of interest
#' @examples
#' \dontrun{
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
#'
#' scoxph_fit = scoxph(x, w, y, D, times)
#' scoxph_cate = predict(scoxph_fit, x, times)
#' }
#' @return a scoxph object
#' @export
scoxph = function(x, w, y, D, times){

  input = sanitize_input(x,w,y,D)
  x = input$x
  w = input$w
  y = input$y
  D = input$D

  x_tilde = data.frame(as.numeric(w - 0.5) * cbind(1, x), x)
  x_pred1 = data.frame(0.5 * cbind(1, x), x)
  x_pred0 = data.frame(-0.5 * cbind(1, x), x)

  colnames(x_tilde) <- colnames(x_pred1) <- colnames(x_pred0) <- paste0("v", 1:dim(x_tilde)[2])
  formula <- as.formula(paste0("survival::Surv(y, D) ~ ", paste(colnames(x_tilde), sep=" ", collapse = "+")))
  tmpdat <- data.frame(y, D, x_tilde)

  s_fit <- survival::coxph(formula, data = tmpdat)
  bh_dat <- survival::basehaz(s_fit, centered = FALSE)
  index <- findInterval(times, bh_dat$time)
  bh <- bh_dat[index, 1]

  link1 <- exp(as.matrix(x_pred1) %*% s_fit$coefficients)
  link0 <- exp(as.matrix(x_pred0) %*% s_fit$coefficients)

  est_S1_cvd <- exp(-bh)^link1
  est_S0_cvd <- exp(-bh)^link0

  tau_hat <- est_S1_cvd - est_S0_cvd

  ret = list(s_fit = s_fit,
             bh = bh,
             s_beta = s_fit$coefficients,
             tau_hat = tau_hat)

  class(ret) <- "scoxph"
  ret
}

#' predict for scoxph
#'
#' get estimated tau(x) using the trained scoxph model
#'
#' @param object A scoxph object
#' @param newx Covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param times The prediction time of interest
#' @param ... Additional arguments (currently not used)
#'
#' @examples
#' \dontrun{
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
#'
#' scoxph_fit = scoxph(x, w, y, D, times)
#' scoxph_cate = predict(scoxph_fit, x, times)
#' }
#' @return vector of estimated conditional average treatment effects
#' @export
predict.scoxph <- function(object,
                           newx = NULL,
                           times,
                           ...) {
  if (!is.null(newx)) {
    newx <- sanitize_x(newx)

    bh_dat <- survival::basehaz(object$s_fit, centered = FALSE)
    index <- findInterval(times, bh_dat$time)
    bh <- bh_dat[index, 1]

    x_pred1 <- data.frame(0.5, 0.5 * newx, newx)
    x_pred0 <- data.frame(-0.5, -0.5 * newx, newx)

    link1 <- exp(as.matrix(x_pred1) %*% object$s_beta)
    link0 <- exp(as.matrix(x_pred0) %*% object$s_beta)

    est_S1_cvd <- exp(-bh)^link1
    est_S0_cvd <- exp(-bh)^link0

    tau_hat <- est_S1_cvd - est_S0_cvd
  }
  else {
    tau_hat <- object$tau_hat
  }
  return(tau_hat)
}
