# *** Comparison methods ***

# S-learner
cate_sl_coxph <- function(data, data.test, t0){
  fit <- surv_sl_coxph(data$X, data$Y, data$W, data$D, t0 = t0)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_sl_lasso <- function(data, data.test, t0){
  fit <- surv_sl_lasso(data$X, data$Y, data$W, data$D, t0 = t0)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_sl_grf <- function(data, data.test, t0, alpha = 0.05){
  fit <- surv_sl_grf(data$X, data$Y, data$W, data$D, t0 = t0, new.args.grf.nuisance = list(alpha = alpha))
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

# T-learner
cate_tl_coxph <- function(data, data.test, t0){
  fit <- surv_tl_coxph(data$X, data$Y, data$W, data$D, t0 = t0)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_tl_lasso <- function(data, data.test, t0){
  fit <- surv_tl_lasso(data$X, data$Y, data$W, data$D, t0 = t0)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_tl_grf <- function(data, data.test, t0, alpha = 0.05){
  fit <- surv_tl_grf(data$X, data$Y, data$W, data$D, t0 = t0, new.args.grf.nuisance = list(alpha = alpha))
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

# X-learner
cate_xl_lasso <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "Kaplan-Meier"){
  fit <- surv_xl_lasso(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_xl_lasso_sf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "survival.forest"){
  fit <- surv_xl_lasso(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_xl_grf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "Kaplan-Meier", alpha = 0.05){
  fit <- surv_xl_grf(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit,
                     new.args.grf.nuisance = list(alpha = alpha))
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_xl_grf_sf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "survival.forest", alpha = 0.05){
  fit <- surv_xl_grf(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit,
                     new.args.grf.nuisance = list(alpha = alpha))
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_xl_grf_lasso <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "Kaplan-Meier", alpha = 0.05){
  fit <- surv_xl_grf_lasso(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit,
                           new.args.grf.nuisance = list(alpha = alpha))
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_xl_grf_lasso_sf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "survival.forest", alpha = 0.05){
  fit <- surv_xl_grf_lasso(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit,
                           new.args.grf.nuisance = list(alpha = alpha))
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

# R-learner
cate_rl_lasso <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "Kaplan-Meier"){
  fit <- surv_rl_lasso(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_rl_lasso_sf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "survival.forest"){
  fit <- surv_rl_lasso(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_rl_grf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "Kaplan-Meier", alpha = 0.05){
  fit <- surv_rl_grf(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit,
                     new.args.grf.nuisance = list(alpha = alpha))
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_rl_grf_sf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "survival.forest", alpha = 0.05){
  fit <- surv_rl_grf(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit,
                     new.args.grf.nuisance = list(alpha = alpha))
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_rl_grf_lasso <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "Kaplan-Meier", alpha = 0.05){
  fit <- surv_rl_grf_lasso(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit,
                           new.args.grf.nuisance = list(alpha = alpha))
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_rl_grf_lasso_sf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "survival.forest", alpha = 0.05){
  fit <- surv_rl_grf_lasso(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit,
                           new.args.grf.nuisance = list(alpha = alpha))
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

# F-learner
cate_fl_lasso <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "Kaplan-Meier"){
  fit <- surv_fl_lasso(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_fl_lasso_sf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "survival.forest"){
  fit <- surv_fl_lasso(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_fl_grf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "Kaplan-Meier"){
  fit <- surv_fl_grf(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

cate_fl_grf_sf <- function(data, data.test, t0, W.hat = 0.5, cen.fit = "survival.forest"){
  fit <- surv_fl_grf(data$X, data$Y, data$W, data$D, t0 = t0, W.hat = W.hat, cen.fit = cen.fit)
  cate <- as.numeric(unlist(predict(fit, data.test$X)))
  cate
}

# Causal survival forest
cate_csf_probs <- function(data, data.test, t0, alpha = 0.05, W.hat = 0.5, cen.fit = NULL) {
  fit <- grf::causal_survival_forest(X = data$X,
                                     Y = data$Y,
                                     W = data$W,
                                     D = data$D,
                                     W.hat = W.hat,
                                     alpha = alpha,
                                     target = "survival.probability",
                                     horizon = t0)
  cate <- as.numeric(unlist(predict(fit, data.test$X)$predictions))
  cate
}
