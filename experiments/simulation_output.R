# Tabulate simulation results
# Run `simulation.R` to produce `simulation.csv.gz`.
# Table 1, Table 2, appendix Figure 5 and Figure 6 are produced below.
rm(list = ls())
library(xtable); library(ggrepel); library(ggplot2)
formatT <- function(x,y){
  paste0(round(x, 2), " (", round(y, 2), ")")
}
setwd("/Users/yizhe/Desktop/Crystal Xu/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/HTE_experiments/sim_output")
old.names <- c("estimate_coxph_sl",
               "estimate_coxph_tl",
               "estimate_csf_probs",
               "estimate_ipcw_las_grf_xl",
               "estimate_ipcw_las_grf_rl",
               "estimate_lasso_sl",
               "estimate_lasso_tl",
               "estimate_ipcw_lasso_fl",
               "estimate_ipcw_lasso_xl",
               "estimate_ipcw_lasso_rl",
               "estimate_grf_sl",
               "estimate_grf_tl",
               "estimate_ipcw_grf_fl",
               "estimate_ipcw_grf_xl",
               "estimate_ipcw_grf_rl")
new.names <- c("sl_coxph","tl_coxph","csf", "xl_las_grf","rl_las_grf",
               "sl_lasso","tl_lasso","fl_lasso", "xl_lasso","rl_lasso",
               "sl_grf","tl_grf","fl_grf", "xl_grf","rl_grf")

grid = expand.grid(n = 5000,
                   p = 25,
                   n.test = 5000,
                   dgp = c("fcomplex"),
                   p_b = c(1, 25, 25),
                   f_b = c("L", "NL", "NL"),
                   pi = c(0.5),
                   gamma = 1,
                   rho = c(2),
                   cen_scale = c(4),
                   cenM = c("indX"),
                   times = (0.2),
                   stringsAsFactors = FALSE)
grid$f_i <- c(rep("L", 6), rep("NL", 3))
grid$p_i <- rep(c(1, 1, 25), 3)
grid <- rbind(grid, grid[2,], grid[5,], grid[8,])
grid[10:12, ]$pi <- c(0.01)               # unbalanced design
grid <- rbind(grid, grid[1,], grid[1,])
grid[13:14, ]$cen_scale <- c(8, 8)        # vary censoring rate (under indX): 30% (default), 70% (early censor), 65%
grid[13, ]$rho <- 1
grid <- rbind(grid, grid[1,])             # vary censoring generating model (= dX)
grid[15, ]$cenM <- "dX"
grid <- rbind(grid, grid[2,], grid[5,], grid[8,], grid[2,], grid[5,], grid[8,])  # vary heterogeneity: sd(CATE)/sd(mu0sp) 0.17, 0.55, 0.9 (baseline)
grid[16:21, ]$gamma <- c(rep(0.46,3), rep(0,3))
grid <- rbind(grid, grid[2,], grid[5,], grid[8,], grid[2,], grid[5,], grid[8,])  # vary event rate
grid[22:27, ]$times <- c(rep(0.02,3), rep(0.001,3))
rownames(grid) <- 1:dim(grid)[1]

df <- NULL
for (k in 1:27){
  dftmp <- NULL
  for (q in 1:5){
    tmp <- read.csv(paste0("040322/ML_grid_",k,"_KMfit_p",q,".csv.gz"), stringsAsFactors = FALSE)
    num <- length(unique(tmp$sim))
    tmp$sim <- seq(num*(q-1)+1, num*q, 1)
    tmp$grid <- k
    dftmp <- rbind(dftmp, tmp)
  }
  df <- rbind(df, dftmp)
}
dfSF <- NULL
for (k in 13:15){
  dftmp <- NULL
  for (q in 1:5){
    tmp <- read.csv(paste0("040322/ML_grid_",k,"_SFfit_p",q,".csv.gz"), stringsAsFactors = FALSE)
    num <- length(unique(tmp$sim))
    tmp$sim <- seq(num*(q-1)+1, num*q, 1)
    tmp$grid <- k
    dftmp <- rbind(dftmp, tmp)
  }
  dfSF <- rbind(dfSF, dftmp)
}
dfSF$cenmodel <- rep("SF", nrow(dfSF))
df$cenmodel <- rep("KM", nrow(df))
df <- rbind(df, dfSF)

gridx <- rbind(grid, grid[13:15,])
grid.df <- data.frame(grid = 1:30, gridx[ ,c(5:dim(gridx)[2])])
df <- merge(x = df, y = grid.df, by = "grid", all.x=T); dim(df)
df <- df[, !colnames(df) %in% c("grid")]
#df$sim <- sort(rep(1:100, 25*6))
#apply(df[c("pi", "dgp", "p_b", "p_i", "f_b", "f_i")], 2, unique)
df$estimator = new.names[match(df$estimator.name, old.names)]
df$estimator.name = NULL
df$classif.error = 1 - df$classif.rate
df$mse = df$mse * 10
df[df$pi == 0.01, ]$dgp <- "unbalanced"
write.csv(df, "full_results.csv", row.names = F)

# -------------------------------------------------------- Tables & Figures --------------------------------------------------------- #
df$p_b_i <- ifelse(df$p_b == 1 & df$p_i == 1, "1_1",
                       ifelse(df$p_b == 25 & df$p_i == 1, "25_1",
                              ifelse(df$p_b == 25 & df$p_i == 25, "25_25", NA)))
df$f_b_i <- ifelse(df$f_b == "L" & df$f_i == "L", "L_L",
                     ifelse(df$f_b == "NL" & df$f_i == "L", "NL_L",
                              ifelse(df$f_b == "NL" & df$f_i == "NL", "NL_NL", NA)))
np_tmp = aggregate(list(MSE = df$mse,
                           Bias = df$bias,
                           Corr = df$rcorr,
                           Calib = df$calib_coef,
                           Concord = df$concordance,
                           CLF = df$classif.error),
                      by = list(estimator = df$estimator,
                                pi = df$pi,
                                dgp = df$dgp,
                                p_b_i = df$p_b_i,
                                f_b_i = df$f_b_i,
                                gamma = df$gamma,
                                rho = df$rho,
                                cen_scale = df$cen_scale,
                                cenM = df$cenM,
                                times = df$times,
                                cenmodel = df$cenmodel),
                      FUN = median)

# figure
library(cowplot)
newnames <- levels(df$estimator)
np_tmp_long <- data.frame(rbind(as.matrix(np_tmp[,c(1:11, 12)]), as.matrix(np_tmp[,c(1:11, 14)])))
names(np_tmp_long) <- c(names(np_tmp)[1:11], "value")
np_tmp_long$Metric <- c(rep("MSE", dim(np_tmp)[1]), rep("Corr", dim(np_tmp)[1]))
np_tmp_long$value <- as.numeric(np_tmp_long$value)
np_tmp_long$rho <- as.numeric(np_tmp_long$rho)
np_tmp_long$cen_scale <- as.numeric(np_tmp_long$cen_scale)
np_tmp_long$gamma <- as.numeric(np_tmp_long$gamma)
np_tmp_long$times <- as.numeric(np_tmp_long$times)

# --- Unbalanced design
np_tmp_long_unbalanced <- np_tmp_long[np_tmp_long$f_b_i == "L_L" &
                                        (np_tmp_long$dgp == "unbalanced" |
                                           (np_tmp_long$p_b_i == "25_1" & np_tmp_long$times == 0.2 & np_tmp_long$rho == 2 & np_tmp_long$gamma == 1 & np_tmp_long$cen_scale == 4)), ]; dim(np_tmp_long_unbalanced)
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_lasso <- np_tmp_long_unbalanced[np_tmp_long_unbalanced$estimator %in% c("sl_lasso", "tl_lasso", "fl_lasso", "xl_lasso", "rl_lasso"), ]
p1 <- ggplot(np_tmp_long_lasso, aes(x = pi, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Treatment Propensity")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p1

np_tmp_long_unbalanced <- np_tmp_long[np_tmp_long$f_b_i == "NL_NL" &
                                        (np_tmp_long$dgp == "unbalanced" |
                                           (np_tmp_long$p_b_i == "25_1" & np_tmp_long$times == 0.2 & np_tmp_long$rho == 2 & np_tmp_long$gamma == 1 & np_tmp_long$cen_scale == 4)), ]; dim(np_tmp_long_unbalanced)
np_tmp_long_grf <- np_tmp_long_unbalanced[np_tmp_long_unbalanced$estimator %in% c("sl_grf", "tl_grf", "fl_grf", "xl_grf", "rl_grf"), ]
p2 <- ggplot(np_tmp_long_grf, aes(x = pi, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Treatment Propensity")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p2

np_tmp_long_unbalanced <- np_tmp_long[np_tmp_long$f_b_i == "NL_L" &
                                        (np_tmp_long$dgp == "unbalanced" |
                                           (np_tmp_long$p_b_i == "25_1" & np_tmp_long$times == 0.2 & np_tmp_long$rho == 2 & np_tmp_long$gamma == 1 & np_tmp_long$cen_scale == 4)), ]; dim(np_tmp_long_unbalanced)
np_tmp_long_lasgrf <- np_tmp_long_unbalanced[np_tmp_long_unbalanced$estimator %in% c("xl_las_grf", "rl_las_grf"), ]
p3 <- ggplot(np_tmp_long_lasgrf, aes(x = pi, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Treatment Propensity")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p3

png("Performance_by_unbalanced.png", width = 11, height = 10, units = 'in', res = 800)
plot_grid(p1, p2, p3, nrow = 2, align = "v")
dev.off()

# --- Censoring rate
np_tmp_long_CR <- np_tmp_long[np_tmp_long$cenM == "indX" & (np_tmp_long$cen_scale != 4 | np_tmp_long$rho != 2| (np_tmp_long$cen_scale == 4 &
                              np_tmp_long$p_b_i == "1_1" & np_tmp_long$f_b_i == "L_L" & np_tmp_long$dgp == "fcomplex" &
                                np_tmp_long$times == 0.2 & np_tmp_long$gamma == 1)), ]; dim(np_tmp_long_CR)
np_tmp_long_CR$groups <- factor(paste0(np_tmp_long_CR$rho, np_tmp_long_CR$cen_scale, np_tmp_long_CR$cenmodel))
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_lasso <- np_tmp_long_CR[np_tmp_long_CR$estimator %in% c("sl_lasso", "tl_lasso", "fl_lasso", "xl_lasso", "rl_lasso"), ]
p1 <- ggplot(np_tmp_long_lasso, aes(x = groups, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Censoring Rate")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~ Metric, scales = "free_y", ncol=2); p1

np_tmp_long_grf <- np_tmp_long_CR[np_tmp_long_CR$estimator %in% c("sl_grf", "tl_grf", "fl_grf", "xl_grf", "rl_grf"), ]
p2 <- ggplot(np_tmp_long_grf, aes(x = groups, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Censoring Rate")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p2

np_tmp_long_lasgrf <- np_tmp_long_CR[np_tmp_long_CR$estimator %in% c("xl_las_grf", "rl_las_grf"), ]
p3 <- ggplot(np_tmp_long_lasgrf, aes(x = groups, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Censoring Rate")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p3

png("Performance_by_CR_indX.png", width = 11, height = 10, units = 'in', res = 800)
plot_grid(p1, p2, p3, nrow = 2, align = "v")
dev.off()


np_tmp_long_CM <- np_tmp_long[np_tmp_long$cenM == "dX" | (np_tmp_long$cen_scale == 4 & np_tmp_long$p_b_i == "1_1" &
                                                            np_tmp_long$f_b_i == "L_L" & np_tmp_long$dgp == "fcomplex"&
                                                            np_tmp_long$times == 0.2 & np_tmp_long$gamma == 1), ]; dim(np_tmp_long_CM)
np_tmp_long_CM$groups <- factor(paste0(np_tmp_long_CM$cenM, np_tmp_long_CM$cenmodel))
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_lasso <- np_tmp_long_CM[np_tmp_long_CM$estimator %in% c("sl_lasso", "tl_lasso", "fl_lasso", "xl_lasso", "rl_lasso"), ]
p1 <- ggplot(np_tmp_long_lasso, aes(x = groups, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Censoring Model")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p1

np_tmp_long_grf <- np_tmp_long_CM[np_tmp_long_CM$estimator %in% c("sl_grf", "tl_grf", "fl_grf", "xl_grf", "rl_grf"), ]
p2 <- ggplot(np_tmp_long_grf, aes(x = groups, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Censoring Rate")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p2

np_tmp_long_lasgrf <- np_tmp_long_CM[np_tmp_long_CM$estimator %in% c("xl_las_grf", "rl_las_grf"), ]
p3 <- ggplot(np_tmp_long_lasgrf, aes(x = groups, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Censoring Rate")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p3

png("Performance_by_CR_dX_NLcen.png", width = 11, height = 10, units = 'in', res = 800)
plot_grid(p1, p2, p3, nrow = 2, align = "v")
dev.off()

# heterogeneity
np_tmp_long_hte_LL <- np_tmp_long[np_tmp_long$f_b_i == "L_L" & (np_tmp_long$gamma != 1 | (np_tmp_long$cen_scale == 4 & np_tmp_long$p_b_i == "25_1" &
                                                            np_tmp_long$dgp == "fcomplex"& np_tmp_long$times == 0.2 & np_tmp_long$gamma == 1)), ]; dim(np_tmp_long_hte_LL)
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_lasso <- np_tmp_long_hte_LL[np_tmp_long_hte_LL$estimator %in% c("sl_lasso", "tl_lasso", "fl_lasso", "xl_lasso", "rl_lasso"), ]
np_tmp_long_lasso$gamma <- as.factor(np_tmp_long_lasso$gamma)
p1 <- ggplot(np_tmp_long_lasso, aes(x = gamma, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("HTE")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p1

np_tmp_long_hte_NLL <- np_tmp_long[np_tmp_long$f_b_i == "NL_L" & (np_tmp_long$gamma != 1 | (np_tmp_long$cen_scale == 4 & np_tmp_long$p_b_i == "25_1" &
                                                                                            np_tmp_long$dgp == "fcomplex"& np_tmp_long$times == 0.2 & np_tmp_long$gamma == 1)), ]; dim(np_tmp_long_hte_NLL)
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_grf <- np_tmp_long_hte_NLL[np_tmp_long_hte_NLL$estimator %in% c("sl_grf", "tl_grf", "fl_grf", "xl_grf", "rl_grf"), ]
np_tmp_long_grf$gamma <- as.factor(np_tmp_long_grf$gamma)
p2 <- ggplot(np_tmp_long_grf, aes(x = gamma, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("HTE")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p2

np_tmp_long_hte_NLNL <- np_tmp_long[np_tmp_long$f_b_i == "NL_NL" & (np_tmp_long$gamma != 1 | (np_tmp_long$cen_scale == 4 & np_tmp_long$p_b_i == "25_1" &
                                                                                              np_tmp_long$dgp == "fcomplex"& np_tmp_long$times == 0.2 & np_tmp_long$gamma == 1)), ]; dim(np_tmp_long_hte_NLNL)
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_lasgrf <- np_tmp_long_hte_NLNL[np_tmp_long_hte_NLNL$estimator %in% c("xl_las_grf", "rl_las_grf"), ]
np_tmp_long_lasgrf$gamma <- as.factor(np_tmp_long_lasgrf$gamma)
p3 <- ggplot(np_tmp_long_lasgrf, aes(x = gamma, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("HTE")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p3

png("Performance_by_hte.png", width = 11, height = 10, units = 'in', res = 800)
plot_grid(p1, p2, p3, nrow = 2, align = "v")
dev.off()

# event rate
np_tmp_long_event_LL <- np_tmp_long[np_tmp_long$f_b_i == "L_L" & (np_tmp_long$times != 0.2 | (np_tmp_long$cen_scale == 4 & np_tmp_long$p_b_i == "25_1" &
                                                                                            np_tmp_long$dgp == "fcomplex"& np_tmp_long$times == 0.2 & np_tmp_long$gamma == 1)), ]; dim(np_tmp_long_event_LL)
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_lasso <- np_tmp_long_event_LL[np_tmp_long_event_LL$estimator %in% c("sl_lasso", "tl_lasso", "fl_lasso", "xl_lasso", "rl_lasso"), ]
np_tmp_long_lasso$times <- as.factor(np_tmp_long_lasso$times)
p1 <- ggplot(np_tmp_long_lasso, aes(x = times, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Event rates")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p1

np_tmp_long_event_NLL <- np_tmp_long[np_tmp_long$f_b_i == "NL_L" & (np_tmp_long$times != 0.2 | (np_tmp_long$cen_scale == 4 & np_tmp_long$p_b_i == "25_1" &
                                                                                              np_tmp_long$dgp == "fcomplex"& np_tmp_long$times == 0.2 & np_tmp_long$gamma == 1)), ]; dim(np_tmp_long_event_NLL)
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_grf <- np_tmp_long_event_NLL[np_tmp_long_event_NLL$estimator %in% c("sl_grf", "tl_grf", "fl_grf", "xl_grf", "rl_grf"), ]
np_tmp_long_grf$times <- as.factor(np_tmp_long_grf$times)
p2 <- ggplot(np_tmp_long_grf, aes(x = times, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Event rates")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p2

np_tmp_long_event_NLNL <- np_tmp_long[np_tmp_long$f_b_i == "NL_NL" & (np_tmp_long$times != 0.2 | (np_tmp_long$cen_scale == 4 & np_tmp_long$p_b_i == "25_1" &
                                                                                                np_tmp_long$dgp == "fcomplex"& np_tmp_long$times == 0.2 & np_tmp_long$gamma == 1)), ]; dim(np_tmp_long_event_NLNL)
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_lasgrf <- np_tmp_long_event_NLNL[np_tmp_long_event_NLNL$estimator %in% c("xl_las_grf", "rl_las_grf"), ]
np_tmp_long_lasgrf$times <- as.factor(np_tmp_long_lasgrf$times)
p3 <- ggplot(np_tmp_long_lasgrf, aes(x = times, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Event rates")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2); p3

png("Performance_by_event_rate.png", width = 11, height = 10, units = 'in', res = 800)
plot_grid(p1, p2, p3, nrow = 2, align = "v")
dev.off()











# --- How sensitive each estimator is to function complexity?
np_tmp_long_fcomplex <- np_tmp_long[np_tmp_long$dgp=="fcomplex" & np_tmp_long$f_b_i=="NL_L",]
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_lasso <- np_tmp_long_fcomplex[np_tmp_long_fcomplex$estimator %in% c("sl_lasso", "tl_lasso", "fl_lasso", "xl_lasso", "rl_lasso"), ]
np_tmp_long_grf <- np_tmp_long_fcomplex[np_tmp_long_fcomplex$estimator %in% c("sl_grf", "tl_grf", "fl_grf", "rl_grf"), ]
np_tmp_long_lasgrf <- np_tmp_long_fcomplex[np_tmp_long_fcomplex$estimator %in% c("sl_coxph", "tl_coxph","csf","xl_las_grf", "rl_las_grf"), ]
p1 <- ggplot(np_tmp_long_lasso, aes(x = p_b_i, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Number of Covariates")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2)
p2 <- ggplot(np_tmp_long_grf, aes(x = p_b_i, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Number of Covariates")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2)
p3 <- ggplot(np_tmp_long_lasgrf, aes(x = p_b_i, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Number of Covariates")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2)

png("Performance_by_fcomplex.png", width = 11, height = 10, units = 'in', res = 800)
plot_grid(p1, p2, p3, nrow = 2, align = "v")
dev.off()





















np_tmp_long_fLNL <- np_tmp_long[np_tmp_long$dgp=="fLNL"|(np_tmp_long$dgp=="fcomplex" & np_tmp_long$p_b_i=="1_1"),]
group.colors <- c("red", "blue", "green","orchid1", "orange")
np_tmp_long_lasso <- np_tmp_long_fLNL[np_tmp_long_fLNL$estimator %in% c("sl_lasso", "tl_lasso", "fl_lasso", "xl_lasso", "rl_lasso"), ]
np_tmp_long_grf <- np_tmp_long_fLNL[np_tmp_long_fLNL$estimator %in% c("sl_grf", "tl_grf", "fl_grf", "xl_grf", "rl_grf"), ]
np_tmp_long_lasgrf <- np_tmp_long_fLNL[np_tmp_long_fLNL$estimator %in% c("sl_coxph", "tl_coxph","csf","xl_las_grf", "rl_las_grf"), ]
p1 <- ggplot(np_tmp_long_lasso, aes(x = f_b_i, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Linear/Nonlinear")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2)
p2 <- ggplot(np_tmp_long_grf, aes(x = f_b_i, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Linear/Nonlinear")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2)
p3 <- ggplot(np_tmp_long_lasgrf, aes(x = f_b_i, y = value)) +
  geom_point(aes(color = estimator), size=2) +
  scale_color_manual(values=group.colors)+
  theme_bw() +
  ylab("") +
  xlab("Linear/Nonlinear")+
  scale_fill_continuous(guide = guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position="bottom")+
  facet_wrap(~Metric, scales = "free_y", ncol=2)
png("Performance_by_fLNL.png", width = 11, height = 10, units = 'in', res = 800)
plot_grid(p1, p2, p3, nrow = 2, align = "v")
dev.off()













# --- Overall performance
overall.table = aggregate(list(MSE = DF$mse,
                          Bias = DF$bias,
                          Corr = DF$rcorr,
                          Calib = DF$calib_coef,
                          Concord = DF$concordance,
                          CLF = DF$classif.error),
                     by = list(estimator = DF$estimator),
                     FUN = median)

overall.90p.table = aggregate(list(MSE.90p = DF$mse,
                                   Bias.90p = DF$bias,
                                   CLF.90p = DF$classif.error),
                              by = list(estimator = DF$estimator),
                              FUN = function(x) quantile(x, prob = c(0.9)))

overall.corr.table = aggregate(list(Corr.10p = DF$rcorr,
                                    Calib.10p = DF$calib_coef,
                                    Concord.10p = DF$concordance),
                          by = list(estimator = DF$estimator),
                          FUN = function(x) quantile(x, prob = c(0.1)))

# Overall performance figure
overall_fig <- overall.table[overall.table$Corr >= 0, ]
ggplot(overall_fig, aes(x = MSE, y = Corr, label = estimator)) +
  geom_point(color = "red", size=1) +
  geom_text_repel(size = 3,
                  box.padding = unit(0.2, "lines")) +
  theme_bw() +
  ylab("Corr") +
  xlab("MSE")
ggsave(paste0("Results_02052022/ml results/overall_performance.png"),
       width = 7, height = 4, dpi = 300)

# Overall performance table
overall.sum <- data.frame(overall.table, overall.90p.table[, 2:dim(overall.90p.table)[2]], overall.corr.table[, 2:dim(overall.corr.table)[2]])
tab1 <- overall.sum
tab1 <- data.frame(tab1[, 1],
                   formatT(tab1$MSE, tab1$MSE.90p), formatT(tab1$Bias, tab1$Bias.90p),
                   formatT(tab1$Corr, tab1$Corr.10p), formatT(tab1$CLF, tab1$CLF.90p))
names(tab1)[2:5] <- c("MSE","Bias","Corr","CLF")
print(xtable(tab1), include.rownames = FALSE)           # generate latex table
write.csv(overall.sum, paste0("Results_02052022/ml results/", "overall_performance.csv"), row.names = FALSE)


# --- How sensitive each estimator is to N to P ratio?
DF2 <- DF[!DF$dgp=="CSF2a",]
np_table = aggregate(list(MSE = DF2$mse,
                          Bias = DF2$bias,
                          Corr = DF2$rcorr,
                          Calib = DF2$calib_coef,
                          Concord = DF2$concordance,
                          CLF = DF2$classif.error),
                     by = list(estimator = DF2$estimator,
                               n = DF2$n,
                               p = DF2$p,
                               np = DF2$np),
                     FUN = median)

# N to P ratio figure
newnames <- levels(DF$estimator)
np_tmp <- np_table[np_table$estimator %in% newnames[c(1, 15, 16, 24)], c(1:5, 7)]
np_tmp_long <- data.frame(rbind(as.matrix(np_tmp[,1:5]), as.matrix(np_tmp[,c(1:4, 6)])))
names(np_tmp_long) <- c(names(np_tmp)[1:4], "value")
np_tmp_long$Metric <- c(rep("MSE", dim(np_tmp)[1]), rep("Corr", dim(np_tmp)[1]))
np_tmp_long$n <- as.factor(np_tmp_long$n)
np_tmp_long$value <- as.numeric(np_tmp_long$value)
np_tmp_long$np <- factor(np_tmp_long$np,
                         levels = c("2000/100","2000/50","5000/100","2000/25","5000/50",
                                    "2000/15","5000/25","5000/15"))
ggplot(np_tmp_long, aes(x = np, y = value)) +
  geom_point(aes(shape = n, color = estimator), size=2) +
  theme_bw() +
  ylab("Performance") +
  xlab("N to P Ratio")+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~Metric, scales = "free_y", ncol=2)
ggsave(paste0("Results_02052022/ml results/N_to_P_ratio.png"),
       width = 7, height = 3, dpi = 300)

# --- How sensitive each estimator is to censoring rate (specified by "horizon")?
DF2 <- DF
cen_table = aggregate(list(MSE = DF2$mse,
                          Bias = DF2$bias,
                          Corr = DF2$rcorr,
                          Calib = DF2$calib_coef,
                          Concord = DF2$concordance,
                          CLF = DF2$classif.error),
                     by = list(estimator = DF2$estimator,
                               n = DF2$n,
                               p = DF2$p,
                               cen = DF2$horizon),
                     FUN = median)

# N to P ratio figure
newnames <- levels(DF$estimator)
np_tmp <- cen_table[cen_table$estimator %in% newnames[c(16, 17, 18, 19)], ]
np_tmp <- cen_table[cen_table$estimator %in% newnames[c(25, 15, 16, 24)], ]
np_tmp_long <- data.frame(rbind(as.matrix(np_tmp[,1:5]), as.matrix(np_tmp[,c(1:4, 7)])))
names(np_tmp_long) <- c(names(np_tmp)[1:4], "value")
np_tmp_long$Metric <- c(rep("MSE", dim(np_tmp)[1]), rep("Corr", dim(np_tmp)[1]))
np_tmp_long$value <- as.numeric(np_tmp_long$value)
np_tmp_long$cen2 <- ifelse(np_tmp_long$cen=="0.20", "7%",
                          ifelse(np_tmp_long$cen=="0.65", "23%",
                                 ifelse(np_tmp_long$cen=="1.22", "44%", "57%")))
np_tmp_long$cen2 <- factor(np_tmp_long$cen2, levels = c("7%", "23%", "44%", "57%"))
ggplot(np_tmp_long, aes(x = cen2, y = value)) +
  geom_point(aes(shape = n, color = estimator), size=2) +
  theme_bw() +
  ylab("Performance") +
  xlab("Censoring rate")+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~Metric, scales = "free_y", ncol=2)
ggsave(paste0("Performance_by_censoring_ER33.png"), width = 7, height = 3, dpi = 300)

# low event rate
np_tmp_long$cen2 <- ifelse(np_tmp_long$cen=="0.02", "16%",
                           ifelse(np_tmp_long$cen=="0.05", "40%",
                                  ifelse(np_tmp_long$cen=="0.08", "61%", "68%")))
np_tmp_long$cen2 <- factor(np_tmp_long$cen2, levels = c("16%", "40%", "61%", "68%"))
ggplot(np_tmp_long, aes(x = cen2, y = value)) +
  geom_point(aes(shape = n, color = estimator), size=2) +
  theme_bw() +
  ylab("Performance") +
  xlab("Censoring rate")+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~Metric, scales = "free_y", ncol=2)
ggsave(paste0("Performance_by_censoring_ER9.png"), width = 7, height = 3, dpi = 300)


# --- How each estimator perform under different dgps?
DF3 <- DF[substr(DF$estimator, 1, 3) != "wcf",]
dgp_table = aggregate(list(MSE = DF3$mse,
                          Bias = DF3$bias,
                          Corr = DF3$rcorr,
                          Calib = DF3$calib_coef,
                          Concord = DF3$concordance,
                          CLF = DF3$classif.error),
                     by = list(estimator = DF3$estimator,
                               dgp = DF3$dgp),
                     FUN = median)

# dgp figure
dgp_table_sub <- dgp_table[dgp_table$Corr>=0.5, ]
dgp_table_sub[dgp_table_sub$dgp == "type2",]$dgp <- "CSF2"
cuts <- aggregate(dgp_table_sub$MSE, by=list(dgp_table_sub$dgp), FUN = median)$x
dgp_table_sub <- dgp_table_sub[dgp_table_sub$dgp == "CSF2" & dgp_table_sub$MSE <= cuts[1]|
                               dgp_table_sub$dgp == "CSF2a" & dgp_table_sub$MSE <= cuts[2]|
                               dgp_table_sub$dgp == "CSF2b" & dgp_table_sub$MSE <= cuts[3]|
                               dgp_table_sub$dgp == "RLsurv1" & dgp_table_sub$MSE <= cuts[4]|
                               dgp_table_sub$dgp == "RLsurv2" & dgp_table_sub$MSE <= cuts[5]|
                               dgp_table_sub$dgp == "RLsurv3" & dgp_table_sub$MSE <= cuts[6], ]
ggplot(dgp_table_sub, aes(x = MSE, y = Corr, label = estimator)) +
  geom_point(color = "red", size=1) +
  geom_text_repel(size = 3,
                  box.padding = unit(0.2, "lines")) +
  theme_bw() +
  ylab("Corr") +
  xlab("MSE")+
  facet_wrap(~dgp, scales = "free", ncol=3)
ggsave(paste0("Results_02052022/ml results/dgp_top_methods.png"),
       width = 8.5, height = 6, dpi = 300)

# dgp tables (total 6)
dgps <- unique(grid$dgp)
for (i in 1:length(dgps)){
  dgpp <- dgps[i]
  index <- row.names(grid[grid$dgp == dgpp,])
  df <- NULL
  for (k in index){
    tmp <- read.csv(paste0("Results_02052022/ml results/metalearner_grid_",k,"_p1.csv.gz"), stringsAsFactors = FALSE)
    df <- rbind(df, tmp)
  }
  df$sim <- sort(rep(1:100, 25 * 8))
  apply(df[c("n", "p", "n.test", "dgp", "estimator.name")], 2, unique)
  df$estimator = new.names[match(df$estimator.name, old.names)]
  df$estimator.name = NULL
  df$classif.error = 1 - df$classif.rate
  df$mse = df$mse * 10

  df.min = aggregate(list(mse.min = df$mse),
                     by = list(n = df$n,
                               p = df$p,
                               n.test = df$n.test,
                               dgp = df$dgp,
                               sim = df$sim),
                     FUN = min)
  DF = merge(df, df.min)
  DF$mse.excess = DF$mse / DF$mse.min

  # Reorder estimator names
  DF$estimator = factor(DF$estimator, levels = new.names)

  # Tables
  DF.table = aggregate(list(MSE = DF$mse,
                            Bias = DF$bias,
                            Corr = DF$rcorr,
                            Calib = DF$calib_coef,
                            Concord = DF$concordance,
                            CLF = DF$classif.error),
                       by = list(estimator = DF$estimator,
                                 dgp = DF$dgp),
                       FUN = median)

  DF.90p.table = aggregate(list(MSE.90p = DF$mse,
                                Bias.90p = DF$bias,
                                CLF.90p = DF$classif.error),
                           by = list(estimator = DF$estimator,
                                     dgp = DF$dgp),
                           FUN = function(x) quantile(x, prob = c(0.9)))

  DF.corr.table = aggregate(list(Corr.10p = DF$rcorr,
                                 Calib.10p = DF$calib_coef,
                                 Concord.10p = DF$concordance),
                            by = list(estimator = DF$estimator,
                                      dgp = DF$dgp),
                            FUN = function(x) quantile(x, prob = c(0.1)))

  DF.sum <- data.frame(DF.table, DF.90p.table[, 3:dim(DF.90p.table)[2]], DF.corr.table[, 3:dim(DF.corr.table)[2]])
  tab1 <- DF.sum[, names(DF.sum) %in% c("estimator", "dgp", "MSE", "MSE.90p", "Bias", "Bias.90p", "Corr", "Corr.10p", "CLF", "CLF.90p")]
  tab1 <- data.frame(tab1[, names(tab1) %in% c("estimator", "dgp")],
                     formatT(tab1$MSE, tab1$MSE.90p), formatT(tab1$Bias, tab1$Bias.90p),
                     formatT(tab1$Corr, tab1$Corr.10p), formatT(tab1$CLF, tab1$CLF.90p))
  names(tab1)<- c("Estimator", "MSE", "Bias", "Corr", "CLF")

  # generate latex table
  print(xtable(tab1), include.rownames = FALSE)
  write.csv(tab1, paste0("Results_02052022/ml results/", DF.sum$dgp[1], ".csv"), row.names = FALSE)
}


# Tables for each N, P, and dgp
ns <- c(2000, 5000, 10000)
for (i in 1:length(dgps)){
  dgpp <- dgps[i]
  for (j in 1:length(ns)){
    nn <- ns[j]
    index <- row.names(grid[grid$dgp == dgpp & grid$n == nn,])
    df <- NULL
    for (k in index){
      tmp <- read.csv(paste0("Results_02052022/ml results/metalearner_grid_",k,"_p1.csv.gz"), stringsAsFactors = FALSE)
      df <- rbind(df, tmp)
    }
    df <- read.csv(paste0("Results_02052022/ml results/metalearner_grid_",k,"_p1.csv.gz"), stringsAsFactors = FALSE)
    df$sim <- sort(rep(1:100, 25 * 4))
    apply(df[c("n", "p", "n.test", "dgp", "estimator.name")], 2, unique)
    df$estimator = new.names[match(df$estimator.name, old.names)]
    df$estimator.name = NULL
    df$classif.error = 1 - df$classif.rate
    df$mse = df$mse * 10

    df.min = aggregate(list(mse.min = df$mse),
                       by = list(n = df$n,
                                 p = df$p,
                                 n.test = df$n.test,
                                 dgp = df$dgp,
                                 sim = df$sim),
                       FUN = min)
    DF = merge(df, df.min)
    DF$mse.excess = DF$mse / DF$mse.min

    # Reorder estimator names
    DF$estimator = factor(DF$estimator, levels = new.names)

    # Tables
    DF.table = aggregate(list(MSE = DF$mse,
                              Bias = DF$bias,
                              Corr = DF$rcorr,
                              Calib = DF$calib_coef,
                              Concord = DF$concordance,
                              CLF = DF$classif.error),
                         by = list(estimator = DF$estimator,
                                   n = DF$n,
                                   p = DF$p,
                                   n.test = DF$n.test,
                                   dgp = DF$dgp),
                         FUN = median)

    DF.90p.table = aggregate(list(MSE.90p = DF$mse,
                                  Bias.90p = DF$bias,
                                  CLF.90p = DF$classif.error),
                             by = list(estimator = DF$estimator,
                                       n = DF$n,
                                       p = DF$p,
                                       n.test = DF$n.test,
                                       dgp = DF$dgp),
                             FUN = function(x) quantile(x, prob = c(0.9)))

    DF.corr.table = aggregate(list(Corr.10p = DF$rcorr,
                                   Calib.10p = DF$calib_coef,
                                   Concord.10p = DF$concordance),
                              by = list(estimator = DF$estimator,
                                        n = DF$n,
                                        p = DF$p,
                                        n.test = DF$n.test,
                                        dgp = DF$dgp),
                              FUN = function(x) quantile(x, prob = c(0.1)))

    DF.sum <- data.frame(DF.table, DF.90p.table[, 6:dim(DF.90p.table)[2]], DF.corr.table[, 6:dim(DF.corr.table)[2]])
    tab1 <- DF.sum[, names(DF.sum) %in% c("estimator","dgp", "n","p", "MSE", "MSE.90p", "Bias", "Bias.90p", "Corr", "Corr.10p", "CLF", "CLF.90p")]
    tab1 <- cbind(tab1[, names(tab1) %in% c("estimator","dgp", "n", "p")],
                  formatT(tab1$MSE, tab1$MSE.90p), formatT(tab1$Bias, tab1$Bias.90p),
                  formatT(tab1$Corr, tab1$Corr.10p), formatT(tab1$CLF, tab1$CLF.90p))
    names(tab1)[5:8] <- c("MSE","Bias","Corr","CLF"); tab1
    #tab1 <- cbind(tab1[tab1$n==2000, !names(tab1) %in% c("n", "p")], tab1[tab1$n==5000, !names(tab1) %in% c("estimator", "n", "p")])

    # generate latex table
    #print(xtable(tab1), include.rownames = FALSE)
    #write.csv(DF.sum, paste0("Results_02052022/ml results/", DF.sum$dgp[1], "_n", DF.sum$n[1], ".csv"), row.names = FALSE)
  }
}








# # *** Performance plots ***
# grid = expand.grid(p = c(25, 50, 100),
#                    n.test = 5000,
#                    dgp = c("RLsurv1", "RLsurv2", "RLsurv3", "RLsurv4",  # R-learner - Nie (2021) - obs
#                            "type1", "type2", "type3", "type4"),         # CSF - Cui (2022) - obs
#                    stringsAsFactors = FALSE)
#
# f <- function(x){
#   paste0(x[,3], "_p", x[,1])
# }
# filenames <- f(grid)
#
# findex <- 16; print(filenames[findex])
# for (j in findex){
#   data <- read.csv(paste0("Results_02052022/ml results/", filenames[j], ".csv"))
#   data$estimator <- factor(data$estimator,
#                            levels = c("coxph_sl", "lasso_sl", "gbm_sl", "grf_sl",
#                                       "coxph_tl", "lasso_tl","gbm_tl","grf_tl",
#                                       "wocf_lasso_xl", "wcf_lasso_xl", "wocf_gbm_xl", "wcf_gbm_xl", "grf_xl",
#                                       "wocf_lasso_fl", "wcf_lasso_fl", "wocf_gbm_fl", "wcf_gbm_fl", "grf_fl",
#                                       "wocf_lasso_rl", "wcf_lasso_rl", "wocf_gbm_rl", "wcf_gbm_rl", "grf_rl",
#                                       "csf_probs", "grf_probs"))
#   data$Learner <- c(rep("X",5), rep("F",5), rep("R",5), rep("S",4), rep("T",4), "CSF", "CF")
#   data$n <- as.factor(data$n)
#   data$MSE <- data$MSE
#
#   # Figure 1 - MSE and Correlation for each n
#   # # data in long format
#   # metrics.col = which(names(data) %in% c("MSE", "Corr"))
#   # data.long = reshape(data,
#   #                         direction = "long",
#   #                         varying = list(names(data)[metrics.col]),
#   #                         v.names = "val",
#   #                         idvar = names(data)[-metrics.col],
#   #                         times = names(data)[metrics.col])
#   # rownames(data.long) = NULL
#   # colnames(data.long)[colnames(data.long) == "time"] = "metric"
#   # ggplot(data.long, aes(x = estimator, y = val, color = Learner, shape=N)) +
#   #   geom_point() +
#   #   theme_bw() +
#   #   ylab("Performance") +
#   #   xlab("Method")+
#   #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   #   facet_wrap(~metric, scales = "free_y", ncol=1)
#   # ggsave(paste0("Results_02052022/ml results/",filenames[j], "_mse_corr.png"), width = 6, height = 5, dpi = 300)
#   datap <- data[data$n==5000, ]
#   ggplot(datap, aes(x = MSE, y = Corr, label = estimator)) +
#     geom_point(color = "red", size=1) +
#     geom_text_repel(size = 3,
#                     box.padding = unit(0.2, "lines")) +
#     theme_bw() +
#     ylab("Corr") +
#     xlab("MSE")
#   ggsave(paste0("Results_02052022/ml results/",filenames[j], "_mse_corr.png"),
#           width = 8, height = 4, dpi = 300)
#
#
#
#
#
#
#
#   # Figure 2 - Bias and CLF for each n (supplement)
#   # data in long format
#   metrics.col = which(names(data) %in% c("Bias", "CLF"))
#   data.long = reshape(data,
#                       direction = "long",
#                       varying = list(names(data)[metrics.col]),
#                       v.names = "val",
#                       idvar = names(data)[-metrics.col],
#                       times = names(data)[metrics.col])
#   rownames(data.long) = NULL
#   colnames(data.long)[colnames(data.long) == "time"] = "metric"
#   names(data.long)[2] <- c("N")
#
#   ggplot(data.long, aes(x = estimator, y = val, color = Learner, shape=N)) +
#     geom_errorbar(aes(ymin=val, ymax=val), width=.1) +
#     geom_point() +
#     theme_bw() +
#     ylab("Performance") +
#     xlab("Method")+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     facet_wrap(~metric, scales = "free_y", ncol=1)
#   ggsave(paste0("Results_02052022/ml results/",filenames[j], "_bias_clf.png"), width = 6, height = 5, dpi = 300)
# }

