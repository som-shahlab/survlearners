
##---------------------- Apply Cox PH method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021     -------------------##
library(survlearners)
library(stringr)
library(ggplot2)
library(ggpubr)
library(grid)
setwd("/Users/yizhe/Desktop/Crystal Xu/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data")

# orginal PCE - added S0 at 5 or 10 years
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

# Combine baseline variables and outcomes into one data frame
blvars <- read.csv("13_baseline_vars.csv")
outcome <- read.csv("outcome.csv")
out_cvd <- outcome[, 1:3]
tmpdat <- merge(blvars, out_cvd, by="MASKID", all=T)
tmpdata <- tmpdat[complete.cases(tmpdat), 2:dim(tmpdat)[2]]; dim(tmpdata) # 9206 19
rownames(tmpdata) <- 1:dim(tmpdata)[1]

# using PCE risk as a predictor
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

sprint_cates_bamrs <- NULL
for (W in 0:1){
  sprint_cates <- NULL
  for (cen.fit in c("Kaplan-Meier", "survival.forest")){
    # create random artificial treatment assignment
    tmpdata_sub <- tmpdata[tmpdata$INTENSIVE == W,]
    set.seed(123)
    tmpdata_sub$INTENSIVE <- rbinom(nrow(tmpdata_sub), 1, 0.5)
    index <- sample(1:nrow(tmpdata_sub), ceiling(nrow(tmpdata_sub)*0.7))
    tmpdata_tr <- tmpdata_sub[index,]
    tmpdata_ts <- tmpdata_sub[-index,]
    sprint_cates_cvd <- list(NA, 3)
    for (vartype in 1:3) {
      if (vartype == 1) {
        tmpdata_train <- list(X = as.matrix(tmpdata_tr[ , 3:15]), Y = tmpdata_tr$t_cvds, W = tmpdata_tr$INTENSIVE, D = tmpdata_tr$cvd)
        tmpdata_test <- list(X = as.matrix(tmpdata_ts[, 3:15]), Y = tmpdata_ts$t_cvds, W = tmpdata_ts$INTENSIVE, D = tmpdata_ts$cvd)
      } else if (vartype == 2) {
        tmpdata_train <- list(X = as.matrix(tmpdata_tr[ , names(tmpdata_tr) %in% c("pcescore", "SUB_CKD")]),
                              Y = tmpdata_tr$t_cvds, W = tmpdata_tr$INTENSIVE, D = tmpdata_tr$cvd)
        tmpdata_test <- list(X = as.matrix(tmpdata_ts[, names(tmpdata_ts) %in% c("pcescore", "SUB_CKD")]),
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

    # ------------------ Statistical Learning Approaches ---------------------#
    # Others
    CSF <- cate_csf_probs(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, alpha = 0.01)
    XGL <- cate_xl_grf_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit, alpha = 0.01)
    RGL <- cate_rl_grf_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit, alpha = 0.01)

    # lasso based
    SL <- cate_sl_lasso(tmpdata_train, tmpdata_test, t0 = t50)
    TL <- cate_tl_lasso(tmpdata_train, tmpdata_test, t0 = t50)
    FL <- cate_fl_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit)
    XLL <- cate_xl_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit)
    RLL <- cate_rl_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit)

    # grf based
    SG <- cate_sl_grf(tmpdata_train, tmpdata_test, t0 = t50, alpha = 0.01)
    TG <- cate_tl_grf(tmpdata_train, tmpdata_test, t0 = t50, alpha = 0.01)
    FG <- cate_fl_grf(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit)
    XGG <- cate_xl_grf(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit, alpha = 0.01)
    RGG <- cate_rl_grf(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit, alpha = 0.01)

    # save results
    sprint_cates_cvd[[vartype]] <- data.frame(rep(W, length(tmpdata_test$Y)), rep(cen.fit, length(tmpdata_test$Y)),
                                              TL, XGG, XGL, RGL, CSF,  SL, FL, XLL, RLL, SG, TG, FG, RGG)
    print(vartype)
    }

    # save temporary tables
    sprint_cates <- rbind(sprint_cates, cbind(c(rep("full", nrow(sprint_cates_cvd[[1]])), rep("pce",nrow(sprint_cates_cvd[[2]])), rep("fixed CKD",nrow(sprint_cates_cvd[[3]]))),
                                              rbind(sprint_cates_cvd[[1]], sprint_cates_cvd[[2]], sprint_cates_cvd[[3]])))
    print(c(W, cen.fit))
  }
  sprint_cates_bamrs <- rbind(sprint_cates_bamrs, sprint_cates)
}
names(sprint_cates_bamrs)[1:3] <- c("vartype", "W", "cen.fit")

# MSE: true CATE is 0 everywhere
sprint_cates_bamrs$group <- paste0(sprint_cates_bamrs$vartype, sprint_cates_bamrs$W, sprint_cates_bamrs$cen.fit)
sprint_out <- list(NA, length(unique(sprint_cates_bamrs$group)))
for (q in 1:length(unique(sprint_cates_bamrs$group))){
  cates_tmp <-  sprint_cates_bamrs[sprint_cates_bamrs$group == unique(sprint_cates_bamrs$group)[q],
                                   !names(sprint_cates_bamrs) %in% c("vartype","W","cen.fit","group")]
  RMSEs <- Bias <- rep(NA, dim(cates_tmp)[2])
  for (k in 1:ncol(cates_tmp)){
    RMSEs[k] <- sqrt(mean((cates_tmp[,k])^2))
    Bias[k] <- mean(abs(cates_tmp[,k]))
  }
  out <- data.frame(Estimator = names(cates_tmp), RMSE = round(RMSEs, 4), Bias = round(Bias, 4))
  sprint_out[[q]] <- out[order(out$RMSE),]
}
names(sprint_out) <- unique(sprint_cates_bamrs$group)
write.csv(sprint_out, "../Analyses Stanford Team/Data Application/Raw Results/sprint_zero_test_pce_out.csv", row.names = F)
write.csv(sprint_cates_bamrs, "../Analyses Stanford Team/Data Application/Raw Results/sprint_zero_pce_cates.csv", row.names = F)

# ------------------------------------------- NULL Analysis on ACCORD ------------------------------------------------ #
accord_blvars <- read.csv("ACCORD_blvars.csv")
accord_blvars <- accord_blvars[, names(accord_blvars) %in% names(blvars)]

# using PCE risk as a predictor
accord_blvars$DM <- rep(0, nrow(accord_blvars))
accord_blvars$RXBP <- ifelse(accord_blvars$N_AGENTS > 0, 1, 0)
pcedat <- data.frame(age = accord_blvars$AGE, totchol = accord_blvars$CHR, hdl = accord_blvars$HDL, sysbp = accord_blvars$SBP,
                     dm = accord_blvars$DM, rxbp = accord_blvars$RXBP, cursmoke = accord_blvars$currentsmoker)
pcedat$grp <- ifelse(accord_blvars$FEMALE == 1 & accord_blvars$RACE_BLACK== 1, 1,
                     ifelse(accord_blvars$FEMALE == 0 & accord_blvars$RACE_BLACK== 1, 3,
                            ifelse(accord_blvars$FEMALE == 0 & accord_blvars$RACE_BLACK== 0, 4,
                                   ifelse(accord_blvars$FEMALE == 1 & accord_blvars$RACE_BLACK== 0, 2, NA))))
pcescore <- rep(NA, nrow(pcedat))
for (z in 1:nrow(pcedat)){
  pcescore[z] <- original.model(pcedat[z,])
}
accord_blvars$pcescore <- pcescore

accord_outcome <- read.csv("ACCORD_outcomes.csv")
accord_outcome <- accord_outcome[, 1:3]
accord_test <- merge(accord_blvars, accord_outcome, by="MASKID",all=T)
accord_test <- accord_test[complete.cases(accord_test),]; dim(accord_test) # 4711 14
accord_test$t_cvds <- ceiling(accord_test$t_cvds)

accord_cates_bamrs <- NULL
for (W in 0:1){
  accord_cates <- NULL
  for (cen.fit in c("Kaplan-Meier", "survival.forest")){
    # create random artificial treatment assignment
    tmpdata_sub <- tmpdata[accord_test$INTENSIVE == W,]
    set.seed(123)
    tmpdata_sub$INTENSIVE <- rbinom(dim(tmpdata_sub)[1], 1, 0.5)
    index <- sample(1:nrow(tmpdata_sub), ceiling(nrow(tmpdata_sub)*0.7))
    tmpdata_tr <- tmpdata_sub[index,]
    tmpdata_ts <- tmpdata_sub[-index,]
    accord_cates_cvd <- list(NA, 3)
    for (vartype in 1:3) {
      if (vartype == 1) {
        tmpdata_train <- list(X = as.matrix(tmpdata_tr[ , 2:14]), Y = tmpdata_tr$t_cvds, W = tmpdata_tr$INTENSIVE, D = tmpdata_tr$cvd)
        tmpdata_test <- list(X = as.matrix(tmpdata_ts[, 2:14]), Y = tmpdata_ts$t_cvds, W = tmpdata_ts$INTENSIVE, D = tmpdata_ts$cvd)
      } else if (vartype == 2) {
        tmpdata_train <- list(X = as.matrix(tmpdata_tr[ , names(tmpdata_tr) %in% c("pcescore", "SUB_CKD")]),
                              Y = tmpdata_tr$t_cvds, W = tmpdata_tr$INTENSIVE, D = tmpdata_tr$cvd)
        tmpdata_test <- list(X = as.matrix(tmpdata_ts[, names(tmpdata_tr) %in% c("pcescore", "SUB_CKD")]),
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

    # ------------------ Statistical Learning Approaches ---------------------#
    # Others
    CSF <- cate_csf_probs(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, alpha = 0.01)
    XGL <- cate_xl_grf_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit, alpha = 0.01)
    RGL <- cate_rl_grf_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit, alpha = 0.01)

    # lasso based
    SL <- cate_sl_lasso(tmpdata_train, tmpdata_test, t0 = t50)
    TL <- cate_tl_lasso(tmpdata_train, tmpdata_test, t0 = t50)
    FL <- cate_fl_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit)
    XLL <- cate_xl_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit)
    RLL <- cate_rl_lasso(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit)

    # grf based
    SG <- cate_sl_grf(tmpdata_train, tmpdata_test, t0 = t50, alpha = 0.01)
    TG <- cate_tl_grf(tmpdata_train, tmpdata_test, t0 = t50, alpha = 0.01)
    FG <- cate_fl_grf(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit)
    XGG <- cate_xl_grf(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit, alpha = 0.01)
    RGG <- cate_rl_grf(tmpdata_train, tmpdata_test, t0 = t50, W.hat = W.hat, cen.fit = cen.fit, alpha = 0.01)

    # save results
    # save results
    accord_cates_cvd[[vartype]] <- data.frame(rep(W, length(tmpdata_test$Y)), rep(cen.fit, length(tmpdata_test$Y)),
                                              TL, XGG, XGL, RGL, CSF, SL, FL, XLL, RLL, SG, TG, FG, RGG)
    }

    # save temporary tables
    accord_cates <- rbind(accord_cates, cbind(c(rep("full", nrow(accord_cates_cvd[[1]])), rep("pce",nrow(accord_cates_cvd[[2]])), rep("fixed CKD",nrow(accord_cates_cvd[[3]]))),
                                              rbind(accord_cates_cvd[[1]], accord_cates_cvd[[2]], accord_cates_cvd[[3]])))
    print(c(W, cen.fit))
  }
  accord_cates_bamrs <- rbind(accord_cates_bamrs, accord_cates)
}
names(accord_cates_bamrs)[1:3] <- c("vartype", "W", "cen.fit")

# MSE: true CATE is 0 everywhere
accord_cates_bamrs$group <- paste0(accord_cates_bamrs$vartype, accord_cates_bamrs$W, accord_cates_bamrs$cen.fit)
accord_out <- list(NA, length(unique(accord_cates_bamrs$group)))
for (q in 1:length(unique(accord_cates_bamrs$group))){
  cates_tmp <-  accord_cates_bamrs[accord_cates_bamrs$group == unique(accord_cates_bamrs$group)[q],
                                   !names(accord_cates_bamrs) %in% c("vartype","W","cen.fit","group")]
  RMSEs <- Bias <- rep(NA, dim(cates_tmp)[2])
  for (k in 1:ncol(cates_tmp)){
    RMSEs[k] <- sqrt(mean((cates_tmp[,k])^2))
    Bias[k] <- mean(abs(cates_tmp[,k]))
  }
  out <- data.frame(Estimator = names(cates_tmp), RMSE = round(RMSEs, 4), Bias = round(Bias, 4))
  accord_out[[q]] <- out[order(out$RMSE),]
}
names(accord_out) <- unique(accord_cates_bamrs$group)
write.csv(accord_out, "../Analyses Stanford Team/Data Application/Raw Results/accord_zero_test_pce_out.csv", row.names = F)
write.csv(accord_cates_bamrs, "../Analyses Stanford Team/Data Application/Raw Results/accord_zero_pce_cates.csv", row.names = F)


# create presentable tables
out_full <- rbind(cbind(sprint_out$full0survival.forest[,1:2], sprint_out$`full0Kaplan-Meier`[,1:2],
                        sprint_out$full1survival.forest[,1:2], sprint_out$`full1Kaplan-Meier`[,1:2]),
                  cbind(accord_out$full0survival.forest[,1:2], accord_out$`full0Kaplan-Meier`[,1:2],
                        accord_out$full1survival.forest[,1:2], accord_out$`full1Kaplan-Meier`[,1:2]))
out_pce <- rbind(cbind(sprint_out$pce0survival.forest[,1:2], sprint_out$`pce0Kaplan-Meier`[,1:2],
                       sprint_out$pce1survival.forest[,1:2], sprint_out$`pce1Kaplan-Meier`[,1:2]),
                 cbind(accord_out$pce0survival.forest[,1:2], accord_out$`pce0Kaplan-Meier`[,1:2],
                       accord_out$pce1survival.forest[,1:2], accord_out$`pce1Kaplan-Meier`[,1:2]))

write.csv(out_full, "../Analyses Stanford Team/Data Application/Raw Results/present_NULL_full.csv", row.names = F)
write.csv(out_pce, "../Analyses Stanford Team/Data Application/Raw Results/present_NULL_pce.csv", row.names = F)

library(xtable)
print(xtable(out_full, digits=c(0, rep(c(0, 4), 4))), include.rownames = FALSE)
print(xtable(out_pce, digits=c(0, rep(c(0, 4), 4))), include.rownames = FALSE)

# # Estimated CATEs against pce linear predictor (SPRINT)
# fig.dat <- data.frame(sprint_cates_cvd[[1]][ ,colnames(sprint_cates_cvd[[1]]) %in% c("XGG", "XGL", "RGL", "CSF", "TL")],
#                       sprint_cates_cvd[[2]][ ,colnames(sprint_cates_cvd[[2]]) %in% c("XGG", "XGL", "RGL", "CSF", "TL")],
#                       sprint_cates_cvd[[3]][ ,colnames(sprint_cates_cvd[[3]]) %in% c("XGG", "XGL", "RGL", "CSF", "TL")],
#                       pcescore = tmpdata_ts$pcescore, CKD = tmpdata_ts$SUB_CKD)
# methodnames <- c("Covariates-based TL CATEs", "Covariates-based XGG CATEs", "Covariates-based XGL CATEs", "Covariates-based RGL CATEs", "Covariates-based CSF CATEs",
#                  "PCE-based TL CATEs", "PCE-based XGG CATEs", "PCE-based XGL CATEs", "PCE-based RGL CATEs", "PCE-based CSF CATEs",
#                  "PCE-based TL CATEs", "PCE-based XGG CATEs", "PCE-based XGL CATEs", "PCE-based RGL CATEs", "PCE-based CSF CATEs")
# figs <- list()
# for (z in 1:(ncol(fig.dat)-2)) {
#   figs[[z]] <- ggplot(fig.dat, aes_string(x = "pcescore", y = names(fig.dat)[z])) +
#     geom_point() +
#     ylab(methodnames[z])
# }
# g <- ggarrange(figs[[1]] + rremove("xlab"), figs[[2]] + rremove("xlab"), figs[[3]] + rremove("xlab"), figs[[4]] + rremove("xlab"),
#                figs[[5]] + rremove("xlab"), figs[[6]] + rremove("xlab"), figs[[7]] + rremove("xlab"), figs[[8]] + rremove("xlab"),
#                figs[[9]] + rremove("xlab"), figs[[10]] + rremove("xlab"), figs[[11]] + rremove("xlab"), figs[[12]] + rremove("xlab"),
#                figs[[13]] + rremove("xlab"), figs[[14]] + rremove("xlab"), figs[[15]]+ rremove("xlab"), ncol=5, nrow=3)
# g <- annotate_figure(g, bottom = textGrob("Ten-year PCE Risk", gp = gpar(cex = 0.9)))
# png(paste0("../Analyses Stanford Team/Data Application/Raw Results/sprint_NULL_PCE_", W, "_", cen.fit, ".png"),
#     width = 10, height = 6, units = 'in', res = 300)
# print(g)
# dev.off()
#
# figs <- list()
# for (z in 1:(ncol(fig.dat)-2)) {
#   figs[[z]] <- ggplot(fig.dat, aes_string(x = "CKD", y = names(fig.dat)[z])) +
#     geom_point() +
#     ylab(methodnames[z])
# }
# g <- ggarrange(figs[[1]] + rremove("xlab"), figs[[2]] + rremove("xlab"), figs[[3]] + rremove("xlab"), figs[[4]] + rremove("xlab"),
#                figs[[5]] + rremove("xlab"), figs[[6]] + rremove("xlab"), figs[[7]] + rremove("xlab"), figs[[8]] + rremove("xlab"),
#                figs[[9]] + rremove("xlab"), figs[[10]] + rremove("xlab"), figs[[11]] + rremove("xlab"), figs[[12]] + rremove("xlab"),
#                figs[[13]] + rremove("xlab"), figs[[14]] + rremove("xlab"), figs[[15]]+ rremove("xlab"), ncol=5, nrow=3)
# g <- annotate_figure(g, bottom = textGrob("CKD", gp = gpar(cex = 0.9)))
# png(paste0("../Analyses Stanford Team/Data Application/Raw Results/sprint_NULL_CKD_", W, "_", cen.fit, ".png"),
#     width = 10, height = 6, units = 'in', res = 300)
# print(g)
# dev.off()

# # Estimated CATEs against pce linear predictor (ACCORD)
# fig.dat <- data.frame(accord_cates_cvd[[1]][ ,colnames(accord_cates_cvd[[1]]) %in% c("XGG", "RGL", "CSF", "TL")],
#                       accord_cates_cvd[[2]][ ,colnames(accord_cates_cvd[[2]]) %in% c("XGG", "RGL", "CSF", "TL")],
#                       pcescore = tmpdata_ts$pcescore, SBP = tmpdata_ts$SBP)
# methodnames <- c("Covariates-based TL CATEs", "Covariates-based XGG CATEs", "Covariates-based RGL CATEs", "Covariates-based CSF CATEs",
#                  "PCE-based TL CATEs", "PCE-based XGG CATEs", "PCE-based RGL CATEs", "PCE-based CSF CATEs")
# figs <- list()
# for (z in 1:(ncol(fig.dat)-2)) {
#   figs[[z]] <- ggplot(fig.dat, aes_string(x = "pcescore", y = names(fig.dat)[z])) +
#     geom_point() +
#     ylab(methodnames[z])
# }
# g <- ggarrange(figs[[1]] + rremove("xlab"), figs[[2]] + rremove("xlab"), figs[[3]] + rremove("xlab"), figs[[4]] + rremove("xlab"),
#                figs[[5]] + rremove("xlab"), figs[[6]] + rremove("xlab"), figs[[7]] + rremove("xlab"), figs[[8]] + rremove("xlab"),
#                ncol=4, nrow=2)
# g <- annotate_figure(g, bottom = textGrob("Ten-year PCE Risk", gp = gpar(cex = 0.9)))
# png(paste0("../Analyses Stanford Team/Data Application/Raw Results/accord_NULL_PCE_", W, "_", cen.fit, ".png"),
#     width = 10, height = 5, units = 'in', res = 300)
# print(g)
# dev.off()
#
# figs <- list()
# for (z in 1:(ncol(fig.dat)-2)) {
#   figs[[z]] <- ggplot(fig.dat, aes_string(x = "SBP", y = names(fig.dat)[z])) +
#     geom_point() +
#     ylab(methodnames[z])
# }
# g <- ggarrange(figs[[1]] + rremove("xlab"), figs[[2]] + rremove("xlab"), figs[[3]] + rremove("xlab"), figs[[4]] + rremove("xlab"),
#                figs[[5]] + rremove("xlab"), figs[[6]] + rremove("xlab"), figs[[7]] + rremove("xlab"), figs[[8]] + rremove("xlab"),
#                ncol=4, nrow=2)
# g <- annotate_figure(g, bottom = textGrob("Systolic Blood Pressure", gp = gpar(cex = 0.9)))
# png(paste0("../Analyses Stanford Team/Data Application/Raw Results/accord_NULL_SBP_", W, "_", cen.fit, ".png"),
#     width = 10, height = 5, units = 'in', res = 300)
# print(g)
# dev.off()




