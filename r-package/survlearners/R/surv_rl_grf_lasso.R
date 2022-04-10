#' @title R-learner of grf and lasso
#'
#' @description  R-learner, implemented via the grf package for nuisance parameter estimation and lasso for target parameter
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
#' @param ps The propensity score
#' @param cen_fit The choice of model fitting for censoring
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
#'
#' cate = surv_rl_grf_lasso(X, W, Y, D, times, ps = 0.5, newX = X)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_rl_grf_lasso <- function(X, W, Y, D, times, ps = NULL, alpha = 0.05, cen_fit = "KM", newX = NULL){
  rlasgrf_fit <- rlasgrf(X = X,
                         W = W,
                         Y = Y,
                         D = D,
                         p_hat = ps,
                         alpha = alpha,
                         times = times,
                         cen_fit = cen_fit)
  rlasgrf_est <- predict(object = rlasgrf_fit, newx = newX)
  as.vector(rlasgrf_est)
}
