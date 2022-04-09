#' @title S-learner of lasso
#'
#' @description  S-learner, implemented via glmnet (lasso)
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
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
#' data <- list(X = X, W = W, Y = Y, D = D)
#' data.test <- list(X = X, W = W, Y = Y, D = D)
#'
#' cate = surv_sl_lasso(data, data.test, times)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_sl_lasso <- function(data, data.test, times){

  slasso_fit <- slasso_surv(x = data$X,
                            w = data$W,
                            y = data$Y,
                            D = data$D,
                            times = times)

  pred_S_lasso <- predict(slasso_fit, newx = data.test$X, times = times)
  pred_S_lasso
}
