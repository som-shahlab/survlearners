#' @include utils.R
#'
#' @title S-learner of Cox PH
#'
#' @description  S-learner, implemented via Cox proportional hazard models
#'
#' @param data The training data set
#' @param data.test The testing data set
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
#' data <- list(X = X, W = W, Y = Y, D = D)
#' data.test <- list(X = X, W = W, Y = Y, D = D)
#'
#' scoxph_surv_cate = estimate_coxph_sl(data, data.test, times)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
estimate_coxph_sl <- function(data, data.test, times){

  scoxph_fit <- scoxph(x = data$X,
                       w = data$W,
                       y = data$Y,
                       D = data$D,
                       times = times)

  pred_S_coxph <- predict(scoxph_fit, newx = data.test$X, times = times)
  pred_S_coxph
}
