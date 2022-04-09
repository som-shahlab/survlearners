rm(list = ls())
library(survival)
library(grf)
library(randomForestSRC)
library(survival)
library(gbm)
library(glmnet)
library(stringr)
ccdf <- function(pred, true){
  concordant <- paircum <- 0
  for (i in 1:(length(pred)-1)){
    for (j in (i+1):length(pred)){
      paircum <- paircum + 1
      if(pred[i] >= pred[j] & true[i] >= true[j] | pred[i] < pred[j] & true[i] < true[j]){
        concordant <- concordant + 1
      }else{
        concordant <- concordant
      }
    }
  }
  return(concordant/paircum)
}
source("./R/sprint_parametric_simulation.R")
source("./R/dgps.R")
source("./R/utils.R")
source("./R/scoxph.R")
source("./R/slasso_surv.R")
source("./R/Flasso.R")
source("./R/Fgrf.R")
source("./R/rlasso.R")
source("./R/rgrf.R")
source("./R/rlasgrf.R")
source("./R/comparison_estimators.R")

# *** Comparison methods ***
estimators = list(estimate_coxph_sl = estimate_coxph_sl,
                  estimate_coxph_tl = estimate_coxph_tl,
                  estimate_csf_probs = estimate_csf_probs,
                  estimate_ipcw_las_grf_xl = estimate_ipcw_las_grf_xl,
                  estimate_ipcw_las_grf_rl = estimate_ipcw_las_grf_rl,

                  estimate_lasso_sl = estimate_lasso_sl,
                  estimate_lasso_tl = estimate_lasso_tl,
                  estimate_ipcw_lasso_fl = estimate_ipcw_lasso_fl,
                  estimate_ipcw_lasso_xl = estimate_ipcw_lasso_xl,
                  estimate_ipcw_lasso_rl = estimate_ipcw_lasso_rl,

                  estimate_grf_sl = estimate_grf_sl,
                  estimate_grf_tl = estimate_grf_tl,
                  estimate_ipcw_grf_fl = estimate_ipcw_grf_fl,
                  estimate_ipcw_grf_xl = estimate_ipcw_grf_xl,
                  estimate_ipcw_grf_rl = estimate_ipcw_grf_rl)

# *** Setup ***
out = list()
n.sim = 20      # change to 35 when "unbalanced"
n.mc = 10000

# Simulations scenarios
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
grid[16:21, ]$gamma <- c(rep(0.46, 3), rep(0, 3))
grid <- rbind(grid, grid[2,], grid[5,], grid[8,], grid[2,], grid[5,], grid[8,])  # vary event rate
grid[22:27, ]$times <- c(rep(0.02,3), rep(0.001,3))
rownames(grid) <- 1:dim(grid)[1]

if(length(args <- commandArgs(T))>0){
  stopifnot(length(args)==1)
  i <- as.integer(args[[1]])
  message("running for grid ", i)
}

n = grid$n[i]
p = grid$p[i]
n.test = grid$n.test[i]
pi = grid$pi[i]
dgp = grid$dgp[i]
p_b = grid$p_b[i]; p_i = grid$p_i[i]
f_b = grid$f_b[i]; f_i = grid$f_i[i]
gamma = grid$gamma[i]
rho = grid$rho[i]
cen_scale = grid$cen_scale[i]
cenM = grid$cenM[i]
times = grid$times[i]
an.error.occured <- rep(NA, n.sim)
for (sim in 1:n.sim) {
  tryCatch( {
  print(paste("sim", sim))
  data = generate_tutorial_survival_data(n = n, p = p, p_b = p_b, p_i = p_i, f_b = f_b, f_i = f_i, pi = pi,
                                         gamma = gamma, rho = rho, cen_scale = cen_scale, cenM = cenM, dgp = dgp,
                                         n.mc = 10, times = times)
  data.test = generate_tutorial_survival_data(n = n, p = p, p_b = p_b, p_i = p_i, f_b = f_b, f_i = f_i, pi = pi,
                                              gamma = gamma, rho = rho, cen_scale = cen_scale, cenM = cenM, dgp = dgp,
                                              n.mc = n.mc, times = times)

  data$Y = pmax(rep(0.001, length(data$Y)), data$Y)
  true.catesp = data.test$catesp
  true.catesp.sign = data.test$catesp.sign

  predictions = matrix(NA, n.test, length(estimators))
  estimator.output = list()
  for (j in 1:length(estimators)) {
    estimator.name = names(estimators)[j]
    print(estimator.name)
    predictions[,j] = as.numeric(unlist(estimators[[estimator.name]](data, data.test, ps = pi, cen_fit = "KM",
                                                                     times = times, meta_learner = TRUE)))
    correct.classification = sign(predictions[,j]) == true.catesp.sign

    # calibration slope
    calib_fit <- lm(predictions[,j] ~ true.catesp)

    dfj = data.frame(
      estimator.name = estimator.name,
      mse = mean((predictions[,j] - true.catesp)^2),
      bias = mean(abs(predictions[,j] - true.catesp)),
      rcorr = cor(predictions[,j], true.catesp),
      #taucorr = cor(predictions[,j], true.catesp, method = "kendall"),  # use Kendall's tau for next round run
      calib_coef = calib_fit$coefficients[2],
      concordance = ccdf(predictions[,j], true.catesp),
      classif.rate = mean(correct.classification, na.rm = TRUE) # NA: to ignore X1 < 0.3 in DGP 4.
    )
    dfj$rcorr[is.na(dfj$rcorr)==TRUE] <- 0  # assign correlation to 0 when CATE = ATE
    estimator.output[[j]] = dfj
  }

  # Scatter plot of pred and true CATEs
  if (sim==1){
    newnames <- str_replace_all(names(estimators), "estimate_", "")
    newnames <- str_replace_all(newnames, "ipcw_", "")
    names(estimators) <- names(estimators)
    png(paste0("grid", i, "cen_fit_KM.png"),
        width = 10, height = 6, units = 'in', res = 300)
    par(mfrow=c(3,5),
        oma = c(4,4,0,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    for (z in 1:length(estimators)){
      plot(predictions[,z], true.catesp, main = newnames[z],
           xlim=c(min(true.catesp), max(true.catesp)),
           ylim=c(min(true.catesp), max(true.catesp)),
           axes = FALSE)
      abline(0, 1, col = "red", lwd=1)
      abline(h=mean(true.catesp), col = "blue", lwd=1)
      levels <- round(seq(min(true.catesp), max(true.catesp), by = (max(true.catesp) - min(true.catesp))/5),1)
      axis(side = 1, at=levels, labels = if (z %in% 11:15) levels else FALSE)
      axis(side = 2, at=levels, labels = if (z %in% c(1, 6, 11)) levels else FALSE)
    }
    title(xlab = "Estimated CATE",
          ylab = "True CATE",
          cex.lab = 1.5,
          outer = TRUE, line = 2)
    dev.off()
  }

  df = do.call(rbind, estimator.output)
  df$n = n
  df$p = p
  df$n.test = n.test
  df$dgp = dgp
  df$horizon = times
  df$sim = sim

  out = c(out, list(df))
  }
  , error = function(e) {an.error.occured[sim] <<- TRUE})
}
print(sum(an.error.occured, na.rm = TRUE))
out.df = do.call(rbind, out)
write.csv(out.df, gzfile(paste0("./ML_grid_",i,"_KMfit_p1.csv.gz")), row.names = FALSE)
