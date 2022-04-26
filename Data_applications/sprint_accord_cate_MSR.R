
##---------------------- Apply Cox PH method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021     -------------------##
library(survlearners)
library(stringr)
library(grf)
library(ggplot2)
library(ggpubr)
source("/Users/yizhe/Desktop/Crystal Xu/Stanford Postdoc/NHLBI R01 Aim 2/risk-vs-hte-R/survival_aipw_score.R")
setwd("/Users/yizhe/Desktop/Crystal Xu/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data")
pe_95_ci <- function(x){
  return(c(mean(x, na.rm = TRUE), quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)))
}
original.link <- function(x) {
  group = x$grp
  age = x$age
  totchol = x$totchol
  hdl = x$hdl
  sysbp = x$sysbp
  dm = x$dm
  rxbp = x$rxbp
  cursmoke = x$cursmoke

  vars = matrix(c(log(age),
                  log(age)^2,
                  log(totchol),
                  log(age)*log(totchol),
                  log(hdl),
                  log(age)*log(hdl),
                  rxbp*log(sysbp),
                  rxbp*log(age)*log(sysbp),
                  (1-rxbp)*log(sysbp),
                  (1-rxbp)*log(age)*log(sysbp),
                  cursmoke,
                  log(age)*cursmoke,
                  dm), ncol=13)

  coefs = list(
    c(17.114, 0, 0.94, 0, -18.920, 4.475, 29.291, -6.432, 27.820, -6.087, 0.691, 0, 0.874),
    c(-29.799, 4.884, 13.54, -3.114, -13.578, 3.149, 2.019, 0, 1.957, 0, 7.574, -1.665, 0.661),
    c(2.469, 0, 0.302, 0, -0.307, 0, 1.916, 0, 1.809, 0, 0.549, 0, 0.645),
    c(12.344, 0, 11.853, -2.664, -7.990, 1.769, 1.797, 0, 1.764, 0, 7.837, -1.795, 0.658))
  mean.risk = c(86.61, -29.18, 19.54, 61.18)
  individual.risk = unlist(vars %*% coefs[[group]])
  link = individual.risk - mean.risk[group]
  return(link)
}
original.model <- function(x, time=10) {
  link = original.link(x)
  relative.risk = exp(link)

  if (time == 10) {
    baseline.survival = c(0.9533, 0.9665, 0.8954, 0.9144)
  }else if (time == 5){
    baseline.survival = c(0.98194, 0.98898, 0.95726, 0.96254)
  }
  group = x$grp
  risk = 1 - baseline.survival[group] ^ relative.risk
  return(risk)
}

# Expected calibration error for treatment heterogeneity (ECETH)
eceth <- function(cate.est, cate.score, g = 10){
  if (length(unique(cate.est)) == 1) {
    warning("CATE predictions are identical, expected squared error of average treatment effect is computed instead")
    error <- mean((cate.score - cate.est)^2)
  } else {
    groups <- cut(cate.est, breaks = unique(quantile(cate.est, probs=seq(0,1,1/g))), include.lowest=TRUE)
    gamma_Delta_hat  <- aggregate(cate.score, by = list(groups), FUN = "mean")$x
    if (length(gamma_Delta_hat) < g){
      warning(paste0(length(gamma_Delta_hat), " groups are used instead of ", g, " due to not enough distinct values"))
    }
    N <- length(cate.score)

    if (N <= g) {
     stop("The number of groups must be smaller than number of predictions")
    }

    gamma_Delta_hat_i <- rep(NA, N)
    for (k in 1:N){
      gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[k]]
      gamma_Delta_hat_i[k] <- (N/g)/(N/g-1)*gamma_hat_i - 1/(N/g-1)*cate.score[k]
    }
    error <- mean((cate.score - cate.est) * (gamma_Delta_hat_i - cate.est))
  }
  return(error)
}

# Combine baseline variables and outcomes into one data frame
blvars <- read.csv("13_baseline_vars.csv")
outcome <- read.csv("outcome.csv")
out_cvd <- outcome[, 1:3]
tmpdat <- merge(blvars, out_cvd, by="MASKID", all=T)
tmpdata <- tmpdat[complete.cases(tmpdat), 2:dim(tmpdat)[2]]; dim(tmpdata) # 9206 13
rownames(tmpdata) <- 1:dim(tmpdata)[1]

# Create pce linear predictor
tmpdata$DM <- rep(0, nrow(tmpdata))
tmpdata$RXBP <- ifelse(tmpdata$N_AGENTS > 0, 1, 0)
pcedat <- data.frame(age = tmpdata$AGE, totchol = tmpdata$CHR, hdl = tmpdata$HDL, sysbp = tmpdata$SBP,
                     dm = tmpdata$DM, rxbp = tmpdata$RXBP, cursmoke = tmpdata$currentsmoker)
pcedat$grp <- ifelse(tmpdata$FEMALE == 1 & tmpdata$RACE_BLACK== 1, 1,
                     ifelse(tmpdata$FEMALE == 0 & tmpdata$RACE_BLACK== 1, 3,
                            ifelse(tmpdata$FEMALE == 0 & tmpdata$RACE_BLACK== 0, 4,
                                   ifelse(tmpdata$FEMALE == 1 & tmpdata$RACE_BLACK== 0, 2, NA))))
pcescore <- rep(NA, nrow(pcedat))
for (z in 1:nrow(pcedat)){
  pcescore[z] <- original.model(pcedat[z,])
}
tmpdata$pcescore <- pcescore
t50 <- floor(365.25*3.26)
set.seed(123)
index <- sample(1:nrow(tmpdata), ceiling(nrow(tmpdata)*0.7))
tmpdata_tr <- tmpdata[index,]
tmpdata_ts <- tmpdata[-index,]

# ------------------ Statistical Learning Approaches ---------------------#
cates_sprint <- autocs <- calib_error <- list(NA, 3)
for (vartype in 1:3) {
  if (vartype == 1) {
    # using actual variables
    tmpdata_train <- list(X = as.matrix(tmpdata_tr[, 2:14]), Y = tmpdata_tr$t_cvds, W = tmpdata_tr$INTENSIVE, D = tmpdata_tr$cvd)
    tmpdata_test <- list(X = as.matrix(tmpdata_ts[, 2:14]), Y = tmpdata_ts$t_cvds, W = tmpdata_ts$INTENSIVE, D = tmpdata_ts$cvd)
  } else if (vartype == 2) {
    # using pce linear predictor + SUB_CKD
    tmpdata_train <- list(X = as.matrix(tmpdata_tr[ , names(tmpdata_tr) %in% c("pcescore","SUB_CKD")]),
                          Y = tmpdata_tr$t_cvds, W = tmpdata_tr$INTENSIVE, D = tmpdata_tr$cvd)
    tmpdata_test <- list(X = as.matrix(tmpdata_ts[ , names(tmpdata_ts) %in% c("pcescore","SUB_CKD")]),
                         Y = tmpdata_ts$t_cvds, W = tmpdata_ts$INTENSIVE, D = tmpdata_ts$cvd)
  } else if (vartype == 3) {
    tmpdata_train <- list(X = as.matrix(tmpdata_tr[ , names(tmpdata_tr) %in% c("pcescore", "SUB_CKD")]),
                          Y = tmpdata_tr$t_cvds, W = tmpdata_tr$INTENSIVE, D = tmpdata_tr$cvd)
    tmpdata_tss <- tmpdata_ts
    tmpdata_tss$SUB_CKD <- rep(0, nrow(tmpdata_tss))
    tmpdata_test <- list(X = as.matrix(tmpdata_tss[, names(tmpdata_tss) %in% c("pcescore", "SUB_CKD")]),
                         Y = tmpdata_tss$t_cvds, W = tmpdata_tss$INTENSIVE, D = tmpdata_tss$cvd)
  }
  W.hat <- mean(tmpdata_train$W)

  # T-learner with lasso
  TL <- cate_tl_lasso(tmpdata_train, tmpdata_test, t0 = t50)

  # X-learner of GRF
  XGG <- cate_xl_grf(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = "survival.forest", alpha = 0.01)

  # X-learner of GRF
  XGL <- cate_xl_grf_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = "survival.forest", alpha = 0.01)

  # R-learner of grf and lasso
  RGL <- cate_rl_grf_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = "survival.forest", alpha = 0.01)

  # CSF
  CSF <- cate_csf_probs(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = "survival.forest", alpha = 0.01)

  # save results
  cates_sprint[[vartype]] <- data.frame(TL, XGG, XGL, RGL, CSF)

  # --------------------------------------- Evaluations -------------------------------------------#
  # RATE
  obs_cate_cvd <- causal_survival_forest(X = tmpdata_test$X,
                                         Y = tmpdata_test$Y,
                                         D = tmpdata_test$D,
                                         W = tmpdata_test$W,
                                         W.hat = mean(tmpdata_test$W),
                                         alpha = 0.01,
                                         target="survival.probability",
                                         horizon = t50)

  set.seed(123)
  TL_priority <- cut(TL, breaks = quantile(TL, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  TL_autoc <- rank_average_treatment_effect(obs_cate_cvd, TL_priority, R = 1000)

  XGG_priority <- cut(XGG, breaks = quantile(XGG, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  XGG_autoc <- rank_average_treatment_effect(obs_cate_cvd, XGG_priority, R = 1000)

  XGL_priority <- cut(XGL, breaks = quantile(XGL, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  XGL_autoc <- rank_average_treatment_effect(obs_cate_cvd, XGL_priority, R = 1000)

  if (length(unique(RGL)) == 1) {
    RGL_autoc <- data.frame(estimate = 0, std.err = 0)
  } else {
    set.seed(123)
    RGL_priority <- cut(RGL, breaks = quantile(RGL, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
    RGL_autoc <- rank_average_treatment_effect(obs_cate_cvd, RGL_priority, R = 1000)
  }

  CSF_priority <- cut(CSF, breaks = quantile(CSF, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  CSF_autoc <- rank_average_treatment_effect(obs_cate_cvd, CSF_priority, R = 1000)

  # summary table of AUTOC
  autoc <- data.frame(names(cates_sprint[[vartype]]),
                      rbind(c(TL_autoc$estimate, TL_autoc$estimate-1.96 * TL_autoc$std.err, TL_autoc$estimate+1.96 * TL_autoc$std.err),
                            c(XGG_autoc$estimate, XGG_autoc$estimate-1.96 * XGG_autoc$std.err, XGG_autoc$estimate+1.96 * XGG_autoc$std.err),
                            c(XGL_autoc$estimate, XGL_autoc$estimate-1.96 * XGL_autoc$std.err, XGL_autoc$estimate+1.96 * XGL_autoc$std.err),
                            c(RGL_autoc$estimate, RGL_autoc$estimate-1.96 * RGL_autoc$std.err, RGL_autoc$estimate+1.96 * RGL_autoc$std.err),
                            c(CSF_autoc$estimate, CSF_autoc$estimate-1.96 * CSF_autoc$std.err, CSF_autoc$estimate+1.96 * CSF_autoc$std.err)))
  autocs[[vartype]] <- data.frame(autoc[,1], round(autoc[, 2:4], 4))
  names(autocs[[vartype]]) <- c("Estimator", "AUTOC", "LB", "UB")

  # HTE calibration error
  nboot <- 100
  g <- 10    # num of bins
  calib <- matrix(NA, nboot, 5)
  Gamma <- rep(NA, nboot)
  for (B in 1:nboot){
    idx <- unique(sample(rownames(tmpdata_test$X), dim(tmpdata_test$X)[1], replace = TRUE))
    tmp_test <- list(X = tmpdata_test$X[rownames(tmpdata_test$X) %in% idx,],
                     Y = tmpdata_test$Y[rownames(tmpdata_test$X) %in% idx],
                     D = tmpdata_test$D[rownames(tmpdata_test$X) %in% idx],
                     W = tmpdata_test$W[rownames(tmpdata_test$X) %in% idx])
    unique_times <- sort(unique(tmp_test$Y))

    # survival curves
    S1 <- survival_forest(X = tmp_test$X[tmp_test$W==1,],
                          Y = tmp_test$Y[tmp_test$W==1],
                          D = tmp_test$D[tmp_test$W==1],
                          prediction.type = "Nelson-Aalen",
                          failure.times = unique_times)
    surv1 <- predict(S1, tmp_test$X, failure.times = unique_times)$predictions

    S0 <- survival_forest(X = tmp_test$X[tmp_test$W==0,],
                          Y = tmp_test$Y[tmp_test$W==0],
                          D = tmp_test$D[tmp_test$W==0],
                          prediction.type = "Nelson-Aalen",
                          failure.times = unique_times)
    surv0 <- predict(S0, tmp_test$X, failure.times = unique_times)$predictions

    C1 <- survival_forest(X = tmp_test$X[tmp_test$W==1,],
                          Y = tmp_test$Y[tmp_test$W==1],
                          D = 1 - tmp_test$D[tmp_test$W==1],
                          prediction.type = "Nelson-Aalen",
                          failure.times = unique_times)
    cen1 <- predict(C1, tmp_test$X, failure.times = unique_times)$predictions

    C0 <- survival_forest(X = tmp_test$X[tmp_test$W==0,],
                          Y = tmp_test$Y[tmp_test$W==0],
                          D = 1 - tmp_test$D[tmp_test$W==0],
                          prediction.type = "Nelson-Aalen",
                          failure.times = unique_times)
    cen0 <- predict(C0, tmp_test$X, failure.times = unique_times)$predictions

    # AIPW score
    Gamma_i <- compute_scores(treatment = tmp_test$W,
                              study_time = tmp_test$Y,
                              status = tmp_test$D,
                              propensity_score = rep(mean(tmp_test$W), nrow(tmp_test$X)),
                              treated_survival_curve = surv1,
                              treated_censoring_survival_curve = cen1,
                              control_survival_curve = surv0,
                              control_censoring_survival_curve = cen0,
                              unique_times = unique_times,
                              end_time = t50,
                              rmst = FALSE)
    Gamma[B] <- mean(Gamma_i)

    # TL calibration error
    TL_cate_sub <- TL[rownames(tmpdata_test$X) %in% idx]
    calib[B,1] <- eceth(TL_cate_sub, Gamma_i)

    # XGG calibration error
    XGG_cate_sub <- XGG[rownames(tmpdata_test$X) %in% idx]
    calib[B,2] <- eceth(XGG_cate_sub, Gamma_i)

    # XGL calibration error
    XGL_cate_sub <- XGL[rownames(tmpdata_test$X) %in% idx]
    calib[B,3] <- eceth(XGL_cate_sub, Gamma_i)

    # RGL calibration error
    RGL_cate_sub <- RGL[rownames(tmpdata_test$X) %in% idx]
    calib[B,4] <- eceth(RGL_cate_sub, Gamma_i)

    # CSF calibration error
    CSF_cate_sub <- CSF[rownames(tmpdata_test$X) %in% idx]
    calib[B,5] <- eceth(CSF_cate_sub, Gamma_i)

    print(B)
  }

  # 95% bootstrapped CI
  calib_error[[vartype]] <- data.frame(round(t(apply(calib, 2, FUN = pe_95_ci)), 4))
  names(calib_error[[vartype]]) <- c("error", "LB", "UB")
  rownames(calib_error[[vartype]]) <- c("TL", "XGG", "XGL", "RGL", "CSF")
  calib_error[[vartype]]$error[calib_error[[vartype]]$error < 0] <- 0
  calib_error[[vartype]]$LB[calib_error[[vartype]]$LB < 0] <- 0
}
write.csv(cates_sprint, "../Analyses Stanford Team/Data Application/sprint_cates.csv", row.names = F)
write.csv(autocs, "../Analyses Stanford Team/Data Application/sprint_autocs.csv", row.names = F)
write.csv(calib_error, "../Analyses Stanford Team/Data Application/sprint_caliberror.csv", row.names = F)

# presentation table
sprint_eval_table <- data.frame(rbind(cbind(autocs[[1]][,1], paste0(autocs[[1]]$AUTOC, " (", autocs[[1]]$LB, ", ", autocs[[1]]$UB, ")"),
                                            paste0(calib_error[[1]]$error, " (", calib_error[[1]]$LB, ", ", calib_error[[1]]$UB, ")")),
                                      cbind(autocs[[2]][,1], paste0(autocs[[2]]$AUTOC, " (", autocs[[2]]$LB, ", ", autocs[[2]]$UB, ")"),
                                            paste0(calib_error[[2]]$error, " (", calib_error[[2]]$LB, ", ", calib_error[[2]]$UB, ")")),
                                      cbind(autocs[[3]][,1], paste0(autocs[[3]]$AUTOC, " (", autocs[[3]]$LB, ", ", autocs[[3]]$UB, ")"),
                                            paste0(calib_error[[3]]$error, " (", calib_error[[3]]$LB, ", ", calib_error[[3]]$UB, ")"))))
names(sprint_eval_table) <- c("Method", "AUTOC", "Calibration Error")
write.csv(sprint_eval_table, "../Analyses Stanford Team/Data Application/Raw Results/present_sprint_full.csv", row.names = F)

library(xtable)
print(xtable(sprint_eval_table), include.rownames = FALSE)



# ------------------------------------------- External Validation ACCORD ------------------------------------------------ #
accord_blvars <- read.csv("ACCORD_blvars.csv")
accord_blvars <- accord_blvars[, names(accord_blvars) %in% names(blvars)]
accord_outcome <- read.csv("ACCORD_outcomes.csv")
accord_outcome <- accord_outcome[, 1:3]
accord_test <- merge(accord_blvars, accord_outcome, by="MASKID",all=T)
accord_test <- accord_test[complete.cases(accord_test),]; dim(accord_test) # 4711 14
accord_test$t_cvds <- ceiling(accord_test$t_cvds)

# Estimated CATE on PCE linear predictor
accord_test$DM <- rep(0, nrow(accord_test))
accord_test$RXBP <- ifelse(accord_test$N_AGENTS > 0, 1, 0)
pcedat <- data.frame(age = accord_test$AGE, totchol = accord_test$CHR, hdl = accord_test$HDL, sysbp = accord_test$SBP,
                     dm = accord_test$DM, rxbp = accord_test$RXBP, cursmoke = accord_test$currentsmoker)
pcedat$grp <- ifelse(accord_test$FEMALE == 1 & accord_test$RACE_BLACK== 1, 1,
                     ifelse(accord_test$FEMALE == 0 & accord_test$RACE_BLACK== 1, 3,
                            ifelse(accord_test$FEMALE == 0 & accord_test$RACE_BLACK== 0, 4,
                                   ifelse(accord_test$FEMALE == 1 & accord_test$RACE_BLACK== 0, 2, NA))))
pcescore <- rep(NA, nrow(pcedat))
for (z in 1:nrow(pcedat)){
  pcescore[z] <- original.model(pcedat[z,])
}
accord_test$pcescore <- pcescore

tmpdata_tr <- tmpdata      # entire sprint data
tmpdata_ts <- accord_test

# ------------------ Statistical Learning Approaches ---------------------#
cates_accord <- autocs.ev <- calib_error.ev <- list(NA, 3)
for (vartype in 1:3) {
  if (vartype == 1) {
    # using actual variables
    tmpdata_train <- list(X = as.matrix(tmpdata_tr[, 2:14]), Y = tmpdata_tr$t_cvds, W = tmpdata_tr$INTENSIVE, D = tmpdata_tr$cvd)
    tmpdata_test <- list(X = as.matrix(tmpdata_ts[, 3:15]), Y = tmpdata_ts$t_cvds, W = tmpdata_ts$INTENSIVE, D = tmpdata_ts$cvd)
  } else if (vartype == 2) {
    # using pce linear predictor + SUB_CKD
    tmpdata_train <- list(X = as.matrix(tmpdata_tr[ , names(tmpdata_tr) %in% c("pcescore","SUB_CKD")]),
                          Y = tmpdata_tr$t_cvds, W = tmpdata_tr$INTENSIVE, D = tmpdata_tr$cvd)
    tmpdata_test <- list(X = as.matrix(tmpdata_ts[ , names(tmpdata_ts) %in% c("pcescore","SUB_CKD")]),
                         Y = tmpdata_ts$t_cvds, W = tmpdata_ts$INTENSIVE, D = tmpdata_ts$cvd)
  } else if (vartype == 3) {
    tmpdata_train <- list(X = as.matrix(tmpdata_tr[ , names(tmpdata_tr) %in% c("pcescore", "SUB_CKD")]),
                          Y = tmpdata_tr$t_cvds, W = tmpdata_tr$INTENSIVE, D = tmpdata_tr$cvd)
    tmpdata_tss <- tmpdata_ts
    tmpdata_tss$SUB_CKD <- rep(0, nrow(tmpdata_tss))
    tmpdata_test <- list(X = as.matrix(tmpdata_tss[, names(tmpdata_tss) %in% c("pcescore", "SUB_CKD")]),
                         Y = tmpdata_tss$t_cvds, W = tmpdata_tss$INTENSIVE, D = tmpdata_tss$cvd)
  }
  W.hat <- mean(tmpdata_train$W)

  # T-learner with lasso
  TL <- cate_tl_lasso(tmpdata_train, tmpdata_test, t0 = t50)

  # X-learner of GRF
  XGG <- cate_xl_grf(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = "survival.forest", alpha = 0.01)

  # X-learner of GRL
  XGL <- cate_xl_grf(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = "survival.forest", alpha = 0.01)

  # R-learner of grf and lasso
  RGL <- cate_rl_grf_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = "survival.forest", alpha = 0.01)

  # CSF
  CSF <- cate_csf_probs(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = "survival.forest", alpha = 0.01)

  # save results
  cates_accord[[vartype]] <- data.frame(TL, XGG, XGL, RGL, CSF)

  # --------------------------------------- Evaluations -------------------------------------------#
  # RATE
  obs_cate_cvd <- causal_survival_forest(X = tmpdata_test$X,
                                         Y = tmpdata_test$Y,
                                         D = tmpdata_test$D,
                                         W = tmpdata_test$W,
                                         W.hat = mean(tmpdata_test$W),
                                         alpha = 0.01,
                                         target="survival.probability",
                                         horizon = t50)

  if (length(unique(RGL)) == 1) {
    RGL_autoc <- data.frame(estimate = 0, std.err = 0)
  } else {
    set.seed(123)
    RGL_priority <- cut(RGL, breaks = quantile(RGL, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
    RGL_autoc <- rank_average_treatment_effect(obs_cate_cvd, RGL_priority, R = 1000)
  }

  if (length(unique(TL)) == 1) {
    TL_autoc <- data.frame(estimate = 0, std.err = 0)
  } else {
    set.seed(123)
    TL_priority <- cut(TL, breaks = quantile(TL, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
    TL_autoc <- rank_average_treatment_effect(obs_cate_cvd, TL_priority, R = 1000)
  }

  set.seed(123)
  XGG_priority <- cut(XGG, breaks = quantile(XGG, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  XGG_autoc <- rank_average_treatment_effect(obs_cate_cvd, XGG_priority, R = 1000)

  XGL_priority <- cut(XGL, breaks = quantile(XGL, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  XGL_autoc <- rank_average_treatment_effect(obs_cate_cvd, XGL_priority, R = 1000)

  CSF_priority <- cut(CSF, breaks = quantile(CSF, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  CSF_autoc <- rank_average_treatment_effect(obs_cate_cvd, CSF_priority, R = 1000)

  # summary table of AUTOC
  autoc <- data.frame(names(cates_accord[[vartype]]),
                      rbind(c(TL_autoc$estimate, TL_autoc$estimate-1.96 * TL_autoc$std.err, TL_autoc$estimate+1.96 * TL_autoc$std.err),
                            c(XGG_autoc$estimate, XGG_autoc$estimate-1.96 * XGG_autoc$std.err, XGG_autoc$estimate+1.96 * XGG_autoc$std.err),
                            c(XGL_autoc$estimate, XGL_autoc$estimate-1.96 * XGL_autoc$std.err, XGL_autoc$estimate+1.96 * XGL_autoc$std.err),
                            c(RGL_autoc$estimate, RGL_autoc$estimate-1.96 * RGL_autoc$std.err, RGL_autoc$estimate+1.96 * RGL_autoc$std.err),
                            c(CSF_autoc$estimate, CSF_autoc$estimate-1.96 * CSF_autoc$std.err, CSF_autoc$estimate+1.96 * CSF_autoc$std.err)))
  autocs.ev[[vartype]] <- data.frame(autoc[,1], round(autoc[, 2:4], 4))
  names(autocs.ev[[vartype]]) <- c("Estimator", "AUTOC", "LB", "UB")

  # HTE calibration error
  nboot <- 100
  g <- 10    # num of bins
  calib_accord <- matrix(NA, nboot, 5)
  Gamma <- rep(NA, nboot)
  for (B in 1:nboot){
    idx <- unique(sample(rownames(tmpdata_test$X), nrow(tmpdata_test$X), replace = TRUE))
    tmp_test <- list(X = tmpdata_test$X[rownames(tmpdata_test$X) %in% idx,],
                     Y = tmpdata_test$Y[rownames(tmpdata_test$X) %in% idx],
                     D = tmpdata_test$D[rownames(tmpdata_test$X) %in% idx],
                     W = tmpdata_test$W[rownames(tmpdata_test$X) %in% idx])
    unique_times <- sort(unique(tmp_test$Y))

    # survival curves
    S1 <- survival_forest(X = tmp_test$X[tmp_test$W==1,],
                          Y = tmp_test$Y[tmp_test$W==1],
                          D = tmp_test$D[tmp_test$W==1],
                          prediction.type = "Nelson-Aalen",
                          failure.times = unique_times)
    surv1 <- predict(S1, tmp_test$X, failure.times = unique_times)$predictions

    S0 <- survival_forest(X = tmp_test$X[tmp_test$W==0,],
                          Y = tmp_test$Y[tmp_test$W==0],
                          D = tmp_test$D[tmp_test$W==0],
                          prediction.type = "Nelson-Aalen",
                          failure.times = unique_times)
    surv0 <- predict(S0, tmp_test$X, failure.times = unique_times)$predictions

    C1 <- survival_forest(X = tmp_test$X[tmp_test$W==1,],
                          Y = tmp_test$Y[tmp_test$W==1],
                          D = 1 - tmp_test$D[tmp_test$W==1],
                          prediction.type = "Nelson-Aalen",
                          failure.times = unique_times)
    cen1 <- predict(C1, tmp_test$X, failure.times = unique_times)$predictions

    C0 <- survival_forest(X = tmp_test$X[tmp_test$W==0,],
                          Y = tmp_test$Y[tmp_test$W==0],
                          D = 1 - tmp_test$D[tmp_test$W==0],
                          prediction.type = "Nelson-Aalen",
                          failure.times = unique_times)
    cen0 <- predict(C0, tmp_test$X, failure.times = unique_times)$predictions

    # AIPW score
    Gamma_i <- compute_scores(treatment = tmp_test$W,
                              study_time = tmp_test$Y,
                              status = tmp_test$D,
                              propensity_score = rep(W.hat, dim(tmp_test$X)[1]),
                              treated_survival_curve = surv1,
                              treated_censoring_survival_curve = cen1,
                              control_survival_curve = surv0,
                              control_censoring_survival_curve = cen0,
                              unique_times = unique_times,
                              end_time = t50,
                              rmst = FALSE)
    Gamma[B] <- mean(Gamma_i)

    # TL.ev calibration error
    TL.ev_cate_sub <- TL[rownames(tmpdata_test$X) %in% idx]
    calib_accord[B,1] <- eceth(TL.ev_cate_sub, Gamma_i)

    # XGG.ev calibration error
    XGG.ev_cate_sub <- XGG[rownames(tmpdata_test$X) %in% idx]
    calib_accord[B,2] <- eceth(XGG.ev_cate_sub, Gamma_i)

    # XGL.ev calibration error
    XGL.ev_cate_sub <- XGL[rownames(tmpdata_test$X) %in% idx]
    calib_accord[B,3] <- eceth(XGL.ev_cate_sub, Gamma_i)

    # RGL.ev calibration error
    RGL.ev_cate_sub <- RGL[rownames(tmpdata_test$X) %in% idx]
    calib_accord[B,4] <- eceth(RGL.ev_cate_sub, Gamma_i)

    # CSF.ev calibration error
    CSF.ev_cate_sub <- CSF[rownames(tmpdata_test$X) %in% idx]
    calib_accord[B,5] <- eceth(CSF.ev_cate_sub, Gamma_i)

    print(B)
  }

  # 95% bootstrapped CI
  calib_error.ev[[vartype]] <- data.frame(round(t(apply(calib_accord, 2, FUN = pe_95_ci)),4))
  names(calib_error.ev[[vartype]]) <- c("error", "LB", "UB")
  rownames(calib_error.ev[[vartype]]) <- c("TL", "XGG", "XGL", "RGL", "CSF")
  calib_error.ev[[vartype]]$error[calib_error.ev[[vartype]]$error < 0] <- 0
  calib_error.ev[[vartype]]$LB[calib_error.ev[[vartype]]$LB < 0] <- 0
}
write.csv(cates_accord, "../Analyses Stanford Team/Data Application/Raw Results/accord_cates.csv", row.names = F)
write.csv(autocs.ev, "../Analyses Stanford Team/Data Application/Raw Results/accord_autocs.csv", row.names = F)
write.csv(calib_error.ev, "../Analyses Stanford Team/Data Application/Raw Results/accord_caliberror.csv", row.names = F)

# presentation table
accord_eval_table <- data.frame(rbind(cbind(autocs.ev[[1]][,1], paste0(autocs.ev[[1]]$AUTOC, " (", autocs.ev[[1]]$LB, ", ", autocs.ev[[1]]$UB, ")"),
                                            paste0(calib_error.ev[[1]]$error, " (", calib_error.ev[[1]]$LB, ", ", calib_error.ev[[1]]$UB, ")")),
                                      cbind(autocs.ev[[2]][,1], paste0(autocs.ev[[2]]$AUTOC, " (", autocs.ev[[2]]$LB, ", ", autocs.ev[[2]]$UB, ")"),
                                            paste0(calib_error.ev[[2]]$error, " (", calib_error.ev[[2]]$LB, ", ", calib_error.ev[[2]]$UB, ")")),
                                      cbind(autocs.ev[[3]][,1], paste0(autocs.ev[[3]]$AUTOC, " (", autocs.ev[[3]]$LB, ", ", autocs.ev[[3]]$UB, ")"),
                                            paste0(calib_error.ev[[3]]$error, " (", calib_error.ev[[3]]$LB, ", ", calib_error.ev[[3]]$UB, ")"))))
names(accord_eval_table) <- c("Method", "AUTOC", "Calibration Error")
write.csv(accord_eval_table, "../Analyses Stanford Team/Data Application/Raw Results/present_accord_full.csv", row.names = F)


library(xtable)
print(xtable(accord_eval_table), include.rownames = FALSE)




# # Estimated CATEs against pce linear predictor (SPRINT)
# fig.dat <- data.frame(cates_sprint[[1]], cates_sprint[[2]], pcescore = tmpdata_ts$pcescore, SBP = tmpdata_ts$SBP)
# methodnames <- c("Covariates-based TL CATEs", "Covariates-based XGG CATEs", "Covariates-based RGL CATEs", "Covariates-based CSF CATEs",
#                  "PCE-based TL CATEs","PCE-based XGG CATEs", "PCE-based RGL CATEs", "PCE-based CSF CATEs")
# figs <- list()
# for (z in 1:8) {
#   figs[[z]] <- ggplot(fig.dat, aes_string(x = "pcescore", y = names(fig.dat)[z])) +
#     geom_point() +
#     ylab(methodnames[z])
# }
# g <- ggarrange(figs[[1]] + rremove("xlab"), figs[[2]] + rremove("xlab"), figs[[3]] + rremove("xlab"), figs[[4]] + rremove("xlab"),
#                figs[[5]] + rremove("xlab"), figs[[6]] + rremove("xlab"), figs[[7]] + rremove("xlab"), figs[[8]] + rremove("xlab"),
#                ncol=4, nrow=2)
# g <- annotate_figure(g, bottom = textGrob("Ten-year PCE Risk", gp = gpar(cex = 0.9)))
# png(paste0("../Analyses Stanford Team/Data Application/Raw Results/sprint_CATE_PCE.png"),width = 10, height = 5, units = 'in', res = 300)
# print(g)
# dev.off()
#
# figs <- list()
# for (z in 1:8) {
#   figs[[z]] <- ggplot(fig.dat, aes_string(x = "SBP", y = names(fig.dat)[z])) +
#     geom_point() +
#     ylab(methodnames[z])
# }
# g <- ggarrange(figs[[1]] + rremove("xlab"), figs[[2]] + rremove("xlab"), figs[[3]] + rremove("xlab"), figs[[4]] + rremove("xlab"),
#                figs[[5]] + rremove("xlab"), figs[[6]] + rremove("xlab"), figs[[7]] + rremove("xlab"), figs[[8]] + rremove("xlab"),
#                ncol=4, nrow=2)
# g <- annotate_figure(g, bottom = textGrob("Systolic Blood Pressure", gp = gpar(cex = 0.9)))
# png(paste0("../Analyses Stanford Team/Data Application/Raw Results/sprint_CATE_SBP.png"),width = 10, height = 5, units = 'in', res = 300)
# print(g)
# dev.off()

# # Estimated CATEs against pce linear predictor (ACCORD)
# fig.dat <- data.frame(cates_accord[[1]], cates_accord[[2]], pcescore = tmpdata_ts$pcescore, SBP = tmpdata_ts$SBP)
# methodnames <- c("Covariates-based TL CATEs", "Covariates-based XGG CATEs", "Covariates-based RGL CATEs", "Covariates-based CSF CATEs",
#                  "PCE-based TL CATEs","PCE-based XGG CATEs", "PCE-based RGL CATEs", "PCE-based CSF CATEs")
# figs <- list()
# for (z in 1:8) {
#   figs[[z]] <- ggplot(fig.dat, aes_string(x = "pcescore", y = names(fig.dat)[z])) +
#     geom_point() +
#     ylab(methodnames[z])
# }
# g <- ggarrange(figs[[1]] + rremove("xlab"), figs[[2]] + rremove("xlab"), figs[[3]] + rremove("xlab"), figs[[4]] + rremove("xlab"),
#                figs[[5]] + rremove("xlab"), figs[[6]] + rremove("xlab"), figs[[7]] + rremove("xlab"), figs[[8]] + rremove("xlab"),
#                ncol=4, nrow=2)
# g <- annotate_figure(g, bottom = textGrob("Ten-year PCE Risk", gp = gpar(cex = 0.9)))
# png(paste0("../Analyses Stanford Team/Data Application/Raw Results/accord_CATE_PCE.png"),width = 10, height = 5, units = 'in', res = 300)
# print(g)
# dev.off()
#
# figs <- list()
# for (z in 1:8) {
#   figs[[z]] <- ggplot(fig.dat, aes_string(x = "SBP", y = names(fig.dat)[z])) +
#     geom_point() +
#     ylab(methodnames[z])
# }
# g <- ggarrange(figs[[1]] + rremove("xlab"), figs[[2]] + rremove("xlab"), figs[[3]] + rremove("xlab"), figs[[4]] + rremove("xlab"),
#                figs[[5]] + rremove("xlab"), figs[[6]] + rremove("xlab"), figs[[7]] + rremove("xlab"), figs[[8]] + rremove("xlab"),
#                ncol=4, nrow=2)
# g <- annotate_figure(g, bottom = textGrob("Systolic Blood Pressure", gp = gpar(cex = 0.9)))
# png(paste0("../Analyses Stanford Team/Data Application/Raw Results/accord_CATE_SBP.png"),width = 10, height = 5, units = 'in', res = 300)
# print(g)
# dev.off()

# # combined RATE plots with calibration plots
# rate_figs <- list(XGG.ev_autoc, RGL.ev_autoc, CSF.ev_autoc)
# rate_names <- rep(c("X-learner with GRF", "R-learner with GRF and LASSO", "CSF"),2)
# png("../Analyses Stanford Team/Data Application/accord_rate_plots.png", width = 9, height = 6, units = 'in', res = 300)
# par(mfrow=c(2,3),
#     oma = c(1,1,0,0) + 0.1,
#     mar = c(2,1,1,1) + 0.1)
# for (z in 1:length(rate_names)){
#   if (z <= 3) {
#     plot(rate_figs[[z]], main = rate_names[z])
#   } else if (z==4) {
#     plot(pred_cate_XGG.ev, obs_cate_XGG.ev,
#          xlim = c(min(c(obs_cate_XGG.ev, pred_cate_XGG.ev)), max(c(obs_cate_XGG.ev, pred_cate_XGG.ev))),
#          ylim = c(min(c(obs_cate_XGG.ev, pred_cate_XGG.ev)), max(c(obs_cate_XGG.ev, pred_cate_XGG.ev))))
#     abline(0, 1, col = "red", lwd = 1)
#   } else if (z==5) {
#     plot(pred_cate_RGL.ev, obs_cate_RGL.ev,
#          xlim = c(min(c(obs_cate_RGL.ev, pred_cate_RGL.ev)), max(c(obs_cate_RGL.ev, pred_cate_RGL.ev))),
#          ylim = c(min(c(obs_cate_RGL.ev, pred_cate_RGL.ev)), max(c(obs_cate_RGL.ev, pred_cate_RGL.ev))))
#     abline(0, 1, col = "red", lwd = 1)
#   } else if (z==6) {
#     plot(pred_cate_CSF.ev, obs_cate_CSF.ev,
#          xlim = c(min(c(obs_cate_CSF.ev, pred_cate_CSF.ev)), max(c(obs_cate_CSF.ev, pred_cate_CSF.ev))),
#          ylim = c(min(c(obs_cate_CSF.ev, pred_cate_CSF.ev)), max(c(obs_cate_CSF.ev, pred_cate_CSF.ev))))
#     abline(0, 1, col = "red", lwd = 1)
#   }
# }
# dev.off()

# # combined RATE plots with calibration plots
# rate_figs <- list(XGG_autoc, RGL_autoc, CSF_autoc)
# rate_names <- rep(c("X-learner with GRF", "R-learner with GRF and LASSO", "CSF"),2)
# png("../Analyses Stanford Team/Data Application/sprint_rate_plots.png", width = 9, height = 6, units = 'in', res = 300)
# par(mfrow=c(2,3),
#     oma = c(1,1,0,0) + 0.1,
#     mar = c(2,1,1,1) + 0.1)
# for (z in 1:length(rate_names)){
#   if (z <= 3) {
#     plot(rate_figs[[z]], main = rate_names[z])
#   } else if (z==4) {
#     plot(XGG_data$pred_cate_XGG, XGG_data$obs_cate_XGG,
#          xlim = c(min(c(obs_cate_XGG, pred_cate_XGG)), max(c(obs_cate_XGG, pred_cate_XGG))),
#          ylim = c(min(c(obs_cate_XGG, pred_cate_XGG)), max(c(obs_cate_XGG, pred_cate_XGG))))
#     abline(0, 1, col = "red", lwd = 1)
#   } else if (z==5) {
#     plot(RGL_data$pred_cate_RGL, RGL_data$obs_cate_RGL,
#          xlim = c(min(c(obs_cate_RGL, pred_cate_RGL)), max(c(obs_cate_RGL, pred_cate_RGL))),
#          ylim = c(min(c(obs_cate_RGL, pred_cate_RGL)), max(c(obs_cate_RGL, pred_cate_RGL))))
#     abline(0, 1, col = "red", lwd = 1)
#   } else if (z==6) {
#     plot(CSF_data$pred_cate_CSF, CSF_data$obs_cate_CSF,
#          xlim = c(min(c(obs_cate_CSF, pred_cate_CSF)), max(c(obs_cate_CSF, pred_cate_CSF))),
#          ylim = c(min(c(obs_cate_CSF, pred_cate_CSF)), max(c(obs_cate_CSF, pred_cate_CSF))))
#     abline(0, 1, col = "red", lwd = 1)
#   }
# }
# dev.off()
