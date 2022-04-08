#' @title S-learner of grf
#'
#' @description  S-learner, implemented via survival_forest (grf package)
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
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
#' data <- list(X = X, W = W, Y = Y, D = D)
#' data.test <- list(X = X, W = W, Y = Y, D = D)
#'
#' cate = surv_sl_grf(data, data.test, times)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_sl_grf <- function(data, data.test, times, alpha = 0.05){

  grffit <- survival_forest(cbind(data$W, data$X),
                            traindat$Y,
                            traindat$D,
                            alpha = alpha,
                            prediction.type = "Nelson-Aalen")
  index <- findInterval(times, grffit$failure.times)
  surf1 <- predict(grffit, cbind(rep(1, length(data.test$Y)), data.test$X))$predictions[, index]
  surf0 <- predict(grffit, cbind(rep(0, length(data.test$Y)), data.test$X))$predictions[, index]
  pred_S_grf <- surf1 - surf0
  pred_S_grf
}
