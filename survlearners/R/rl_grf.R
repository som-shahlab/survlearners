# R-learner of GRF
estimate_ipcw_grf_rl <- function(data, data.test, ps = NULL, times=times,
                                 alpha = 0.05, cen_fit = "KM", meta_learner = TRUE){
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  rgrf_fit <- rgrf(x = data$X,
                   w = data$W,
                   y = data$Y,
                   D = data$D,
                   p_hat = ps,
                   alpha = alpha,
                   failure.times = Y.grid,
                   times = times,
                   cen_fit = cen_fit,
                   meta_learner = meta_learner)
  rgrf_est <- predict(object = rgrf_fit, data.test$X, meta_learner = meta_learner)
  rgrf_est
}
