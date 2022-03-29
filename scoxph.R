#' @include utils.R
#'
#' @title S-learner, implemented via glmnet (lasso)
#'
#' @description  S-learner, as proposed by Imai and Ratkovic (2013), implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param D event indicator
#' @param times prediction time of interest
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' slasso_fit = slasso(x, w, y)
#' slasso_est = predict(slasso_fit, x)
#' }
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
  formula <- as.formula(paste0("Surv(y, D) ~ ", paste(colnames(x_tilde), sep=" ", collapse = "+")))
  tmpdat <- data.frame(y, D, x_tilde)
  
  s_fit <- coxph(formula, data = tmpdat)
  bh_dat <- basehaz(s_fit, centered = FALSE)
  index <- findInterval(times, bh_dat$time)
  bh <- bh_dat[index, 1]

  # Adjust the coefficients for treatment and interaction terms by multiply by 2
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
#' get estimated tau(x) using the trained slasso model
#'
#' @param object a slasso object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param ... additional arguments (currently not used)
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' slasso_fit = slasso(x, w, y)
#' slasso_est = predict(slasso_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.scoxph <- function(object,
                           newx = NULL,
                           times,
                           ...) {
  if (!is.null(newx)) {
    newx <- sanitize_x(newx)
    
    bh_dat <- basehaz(object$s_fit, centered = FALSE)
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
