#' @title T-learner of Cox PH
#'
#' @description  T-learner, implemented via Cox proportional hazard models
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
#' cate = surv_tl_coxph(X, W, Y, D, times, newX = X)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_tl_coxph <- function(X, W, Y, D, times, newX = NULL){

  traindat <- data.frame(Y = Y, D = D, W = W, X)
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]

  # Model for W = 1
  coxph_fit1 <- survival::coxph(survival::Surv(Y, D) ~., data = traindat1)
  bh_dat <- survival::basehaz(coxph_fit1, centered = FALSE)
  index <- findInterval(times, bh_dat$time)
  bh <- bh_dat[index, 1]
  est_r1 <- predict(coxph_fit1, newdata = data.frame(newX), type="risk")
  surf1 <- exp(-bh)^est_r1

  # Model for W = 0
  coxph_fit0 <- survival::coxph(survival::Surv(Y, D) ~., data = traindat0)
  bh_dat <- survival::basehaz(coxph_fit0, centered = FALSE)
  index <- findInterval(times, bh_dat$time)
  bh <- bh_dat[index, 1]
  est_r0 <- predict(coxph_fit0, newdata = data.frame(newX), type="risk")
  surf0 <- exp(-bh)^est_r0

  pred_T_coxph <- surf1 - surf0
  pred_T_coxph
}
