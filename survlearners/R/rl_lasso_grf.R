# R-learner of lasso (target) and grf (nuisance)
estimate_ipcw_las_grf_rl <- function(data, data.test, ps = NULL, times=times, k_folds = 10, lambda_choice = "lambda.min",
                                     alpha = 0.05, cen_fit = "KM", meta_learner = TRUE){
  rlasgrf_fit <- rlasgrf(x = data$X,
                         w = data$W,
                         y = data$Y,
                         D = data$D,
                         p_hat = ps,
                         alpha = alpha,\
                         times = times,
                         k_folds = k_folds,
                         lambda_choice = lambda_choice,
                         cen_fit = cen_fit,
                         meta_learner = meta_learner)
  rlasgrf_est <- predict(object = rlasgrf_fit, newx = data.test$X)
  rlasgrf_est
}
