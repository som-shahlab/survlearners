#' @title causal survival forest
#'
#' @description  causal survival forest, implemented via the grf R package
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
#' @param ps The propensity score
#' @param cen_fit The choice of model fitting for censoring
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
#' csf_cate = estimate_csf_probs(data, data.test, times, ps = 0.5)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
estimate_csf_probs <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM") {
  fit <- grf::causal_survival_forest(X = data$X,
                                     Y = data$Y,
                                     W = data$W,
                                     D = data$D,
                                     W.hat = ps,
                                     alpha = alpha,
                                     target="survival.probability",
                                     horizon = times)
  p <- predict(fit, data.test$X)
  p$predictions
}
