#' @title S-learner for survival data implemented via grf::survival_forest
#'
#' @description  S-learner via grf::survival_forest
#'
#' @param data training data
#' @param data.test testing data
#' @param time time at which to predict
#'
#' @export
#'
surv_s_grf <- function(data, data.test, time, alpha = 0.05) {

  m <- grf::survival_forest(cbind(data$W, data$X), data$Y,
                       data$D, alpha = alpha, prediction.type = "Nelson-Aalen",
                       failure.times = seq(min(data$Y), max(data$Y), length.out = 101))

  m1_hat <- predict(m, cbind(1, data.test$X), time)$predictions[,1]
  m0_hat <- predict(m, cbind(0, data.test$X), time)$predictions[,1]
  tau_hat <- m1_hat - m0_hat
}
