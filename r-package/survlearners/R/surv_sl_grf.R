#' @title S-learner of grf
#'
#' @description  S-learner, implemented via survival_forest (grf package)
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
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
#' cate = surv_sl_grf(X, W, Y, D, times, newX = X)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_sl_grf <- function(X, W, Y, D, times, alpha = 0.05, newX = NULL){

  grffit <- grf::survival_forest(cbind(W, X),
                                 Y,
                                 D,
                                 alpha = alpha,
                                 prediction.type = "Nelson-Aalen")
  index <- findInterval(times, grffit$failure.times)
  surf1 <- predict(grffit, cbind(rep(1, length(data.test$Y)), newX))$predictions[, index]
  surf0 <- predict(grffit, cbind(rep(0, length(data.test$Y)), newX))$predictions[, index]
  pred_S_grf <- surf1 - surf0
  pred_S_grf
}
