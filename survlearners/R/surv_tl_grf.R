#' @title T-learner of grf
#'
#' @description  T-learner, implemented via survival_forest (grf package)
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
#' cate = surv_tl_grf(data, data.test, times)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_tl_grf <- function(data, data.test, times, alpha = 0.05){
  # Model for W = 1
  grffit1 <- grf::survival_forest(data$X[dataW==1],
                                  data$Y[dataW==1],
                                  data$D[dataW==1],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
  index <- findInterval(times, grffit1$failure.times)
  surf1 <- predict(grffit1, data.test$X)$predictions[, index]

  # Model for W = 0
  grffit0 <- grf::survival_forest(data$X[dataW==0],
                                  data$Y[dataW==0],
                                  data$D[dataW==0],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
  index <- findInterval(times, grffit0$failure.times)
  surf0 <- predict(grffit0, data.test$X)$predictions[, index]

  pred_T_grf <- surf1 - surf0
  pred_T_grf
}
