#' @title T-learner of Cox PH
#'
#' @description  T-learner, implemented via Cox proportional hazard models
#'
#' @param x The baseline covariates
#' @param w The treatment variable (0 or 1)
#' @param y The follow-up time
#' @param D The event indicator
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
#' cate = surv_tl_coxph(data, data.test, times)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_tl_coxph <- function(data, data.test, times){

  traindat <- data.frame(Y = data$Y, D = data$D, W = data$W, data$X)
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]

  # Model for W = 1
  coxph_fit1 <- coxph(Surv(Y, D) ~., data = traindat1)
  bh_dat <- basehaz(coxph_fit1, centered = FALSE)
  bh <- bh_dat[which.min(abs(bh_dat$time - times)),]$hazard
  est_r1 <- predict(coxph_fit1, newdata = data.frame(data.test$X), type="risk")
  surf1 <- exp(-bh)^est_r1

  # Model for W = 0
  coxph_fit0 <- coxph(Surv(Y, D) ~., data = traindat0)
  bh_dat <- basehaz(coxph_fit0, centered = FALSE)
  bh <- bh_dat[which.min(abs(bh_dat$time - times)),]$hazard
  est_r0 <- predict(coxph_fit0, newdata = data.frame(data.test$X), type="risk")
  surf0 <- exp(-bh)^est_r0

  pred_T_coxph <- surf1 - surf0
  pred_T_coxph
}
