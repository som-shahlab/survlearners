#' @title S-learner of lasso
#'
#' @description  S-learner, implemented via glmnet (lasso)
#'
#' @param X The baseline covariates
#' @param W The treatment variable (0 or 1)
#' @param Y The follow-up time
#' @param D The event indicator
#' @param times The prediction time of interest
#' @param alpha Mix tuning parameter for the elastic net
#' @param k_folds Number of folds for cross validation
#' @param foldid User-supplied foldid. Must have length equal to length(W). If provided, it overrides the k_folds option.
#' @param lambda User-supplied lambda sequence for cross validation
#' @param lambda_choice How to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty_factor User-supplied penalty factor, must be of length the same as number of features in X
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
#'
#' surv_sl_lasso_fit = surv_sl_lasso(X, W, Y, D, times)
#' cate = predict(surv_sl_lasso_fit)
#' }
#' @return a surv_sl_lasso object
#' @export
surv_sl_lasso = function(X, W, Y, D, times,
                         alpha = 1,
                         k_folds = NULL,
                         foldid = NULL,
                         lambda = NULL,
                         lambda_choice = "lambda.min",
                         penalty_factor = NULL){

  input = sanitize_input(X, W, Y, D)
  X = input$X
  W = input$W
  Y = input$Y
  D = input$D

  x_scl = scale(X, center = TRUE, scale = TRUE)
  x_scl = x_scl[,!is.na(colSums(x_scl)), drop = FALSE]

  nobs = nrow(x_scl)
  pobs = ncol(x_scl)

  if (is.null(foldid) || length(foldid) != length(W)) {

    if (!is.null(foldid) && length(foldid) != length(W)) {
      warning("supplied foldid does not have the same length ")
    }

    if (is.null(k_folds)) {
      k_folds = floor(max(3, min(10,length(W)/4)))
    }

    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(k_folds), length = length(W)))
  }

  x_scl_tilde = cbind(as.numeric(2 * W - 1) * cbind(1, x_scl), x_scl)
  x_scl_pred1 = cbind(1, x_scl, x_scl)
  x_scl_pred0 = cbind(0, 0 * x_scl, x_scl)

  if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
    if (!is.null(penalty_factor) && length(penalty_factor) != 2 * pobs + 1) {
      warning("penalty_factor supplied is not 1 plus 2 times the number of columns in X. Using all ones instead.")
    }
    penalty_factor = c(0, rep(1, 2 * pobs))
  }
  x_scl_tilde <- as.matrix(data.frame(x_scl_tilde))
  s_fit <- glmnet::cv.glmnet(x_scl_tilde,
                             survival::Surv(Y, D),
                             family = "cox",
                             foldid = foldid,
                             lambda = lambda,
                             penalty.factor = penalty_factor,
                             standardize = FALSE,
                             alpha = alpha)

  s_beta <- t(as.vector(coef(s_fit, s = lambda_choice)))
  s_beta_adj <- c(0.5 * s_beta[1:(1 + dim(X)[2])], s_beta[(2 + dim(X)[2]):dim(x_scl_tilde)[2]])

  link1 <- exp(x_scl_pred1 %*% s_beta_adj)
  link0 <- exp(x_scl_pred0 %*% s_beta_adj)
  S0_t <- base_surv(fit = s_fit,
                    Y = Y,
                    D = D,
                    X = x_scl_tilde,
                    lambda = s_fit$lambda.min)
  index <- findInterval(times, S0_t$time)
  S0 <- S0_t[index,]$survival
  surv1 <- S0^exp(link1)
  surv0 <- S0^exp(link0)

  tau_hat <- as.numeric(surv1 - surv0)

  ret = list(s_fit = s_fit,
             x_org = x_scl_tilde,
             y_org = Y,
             D_org = D,
             beta_org = s_beta,
             s_beta = s_beta_adj,
             S0_t = S0_t,
             times = times,
             tau_hat = tau_hat,
             lambda_choice = lambda_choice)

  class(ret) <- "surv_sl_lasso"
  ret
}

#' predict for surv_sl_lasso
#'
#' get estimated tau(X) using the trained surv_sl_lasso model
#'
#' @param object A surv_sl_lasso object
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
#'
#' surv_sl_lasso_fit = surv_sl_lasso(X, W, Y, D, times)
#' cate = predict(surv_sl_lasso_fit)
#' }
#'
#' @return vector of estimated conditional average treatment effects
#' @export
predict.surv_sl_lasso <- function(object,
                                  newx = NULL,
                                  times = NULL,
                                  ...) {
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    newx_scl = scale(newx, center = TRUE, scale = TRUE)
    newx_scl = newx_scl[,!is.na(colSums(newx_scl)), drop = FALSE]
    newx_scl_pred1 = cbind(1, newx_scl, newx_scl)
    newx_scl_pred0 = cbind(0, 0 * newx_scl, newx_scl)

    link1 <- exp(newx_scl_pred1 %*% object$s_beta)
    link0 <- exp(newx_scl_pred0 %*% object$s_beta)

    if(is.null(times)){
    times <- object$times
    }
    index <- findInterval(times, object$S0_t$time)
    S0 <- object$S0_t[index,]$survival

    surv1 <- S0_t^exp(link1)
    surv0 <- S0_t^exp(link0)

    tau_hat <- as.numeric(surv1 - surv0)
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
