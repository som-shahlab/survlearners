#' @title S-learner of lasso
#'
#' @description  S-learner, implemented via glmnet (lasso)
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
#'
#' cate = surv_sl_lasso(X, W, Y, D, times, newX = X)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_sl_lasso <- function(X, W, Y, D, times, newX = NULL){

  slasso_fit <- slasso_surv(X = X,
                            W = W,
                            Y = Y,
                            D = D,
                            times = times)

  pred_S_lasso <- predict(slasso_fit, newx = newX, times = times)
  pred_S_lasso
}
