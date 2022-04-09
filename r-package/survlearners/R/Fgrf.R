#' @title F-learner of grf
#'
#' @description  F-learner, implemented via the grf package
#'
#' @param x matrix of covariates
#' @param w vector of treatment indicators (0 or 1)
#' @param y vector of response values
#' @param pscore vector of propensity scores
#' @param num.trees number of trees
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
#' Fgrf_fit = Fgrf(data$X, data$W, data$Y, weight = weight)
#' Fgrf_cate = predict(Fgrf_fit, data.test$X)
#' }
#' @return An Flasso object
#' @export
Fgrf = function(x, w, y, pscore = rep(.5, nrow(x)), num.trees = 2000, weight) {

  # Input sanitization
  x <- as.matrix(x)
  if (nrow(x) != length(w)) {
    stop('nrow(x) does not match length(w)')

  } else if (nrow(x) != length(y)) {
    stop('nrow(x) does not match length(y)')

  } else if (!is.numeric(x)) {
    stop('x must be numeric matrix')

  } else if (!is.numeric(y)) {
    stop('y must be numeric (use 0/1 for binary response)')

  } else if (!is.numeric(w) | length(setdiff(w, 0:1)) > 0) {
    stop('w must be vector of 0s and 1s')

  }

  fit = list(x = x, pscore = pscore, ipcw = weight)
  z = w * y / pscore - (1 - w) * y / (1 - pscore)
  fit$tau_fit <- grf::regression_forest(x, z, sample.weights = weight)
  class(fit) = 'Fgrf'
  fit
}


#' predict for Fgrf
#'
#' get estimated tau(x) using the trained Fgrf model
#'
#' @param object An Fgrf object
#' @param newx Covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
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
#' Fgrf_fit = Fgrf(data$X, data$W, data$Y, weight)
#' Fgrf_cate = predict(Fgrf_fit, data.test$X)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.Fgrf = function(object,
                        newx,
                        ...) {
  return(predict(object$tau_fit, newdata = newx)$predictions)
}
