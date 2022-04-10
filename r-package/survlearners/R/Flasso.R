#' @title F-learner of lasso
#'
#' @description  F-learner, implemented via glmnet (lasso)
#'
#' @param X matrix of covariates
#' @param W vector of treatment indicators (0 or 1)
#' @param Y vector of response values
#' @param pscore vector of propensity scores
#' @param nfolds number of cross-validation
#' @param alpha mixing tuning parameter for lasso
#' @param weight vector of subject level weights
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
#' weight <- rep(1, length(Y))
#' data <- list(X = X, W = W, Y = Y, D = D)
#' data.test <- list(X = X, W = W, Y = Y, D = D)
#'
#' Flasso_fit = Flasso(data$X, data$W, data$Y, weight = weight)
#' Flasso_cate = predict(Flasso_fit, data.test$X)
#' }
#' @return An Flasso object
#' @export
Flasso = function(X, W, Y, weight, pscore = rep(.5, nrow(X)), nfolds = 10, alpha = 1) {

  # Input sanitization
  X = as.matrix(X)
  if (nrow(X) != length(W)) {
    stop('nrow(X) does not match length(W)')

  } else if (nrow(X) != length(Y)) {
    stop('nrow(X) does not match length(Y)')

  } else if (!is.numeric(X)) {
    stop('X must be numeric matrix')

  } else if (!is.numeric(Y)) {
    stop('Y must be numeric (use 0/1 for binary response)')

  } else if (!is.numeric(W) | length(setdiff(W, 0:1)) > 0) {
    stop('W must be vector of 0s and 1s')

  }

  fit = list(X = X, pscore = pscore, ipcw = weight)
  Z = W * Y / pscore - (1 - W) * Y / (1 - pscore)
  fit$tau_fit <- glmnet::cv.glmnet(X, Z, family = "gaussian", weights = weight, nfolds = nfolds, alpha = alpha)
  class(fit) = 'Flasso'
  fit
}


#' predict for Flasso
#'
#' get estimated tau(X) using the trained Flasso model
#'
#' @param object An Flasso object
#' @param newx Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
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
#' weight <- rep(1, length(Y))
#' data <- list(X = X, W = W, Y = Y, D = D)
#' data.test <- list(X = X, W = W, Y = Y, D = D)
#'
#' Flasso_fit = Flasso(data$X, data$W, data$Y, weight)
#' Flasso_cate = predict(Flasso_fit, data.test$X)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.Flasso = function(object,
                          newx,
                          ...) {
  return(as.vector(predict(object$tau_fit, newx = newx, s = "lambda.min")))
}
