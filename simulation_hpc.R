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

# set.seed(123)
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
n.sim = 20
n.mc = 10000

# Simulations scenarios
grid = expand.grid(n = 5000,
                    p = 25,
                    n.test = 5000,
                    p_b = c(1, 25, 25),
                    f_b = c("NL", "L", "NL"),
                    dgp = c("fcomplex", "fLNL"),
                    stringsAsFactors = FALSE)
grid <- grid[c(4:6, 10, 16), ]
grid$p_i <- c(1, 1, 25, rep(1, 2))
grid$f_i <- c(rep("L", 3), "L", "NL")
rownames(grid) <- 1:5

if(length(args <- commandArgs(T))>0){
  stopifnot(length(args)==1)
  i <- as.integer(args[[1]])
  message("running for grid ", i)
}

n = grid$n[i]
p = grid$p[i]
n.test = grid$n.test[i]
dgp = grid$dgp[i]
p_b = grid$p_b[i]; p_i = grid$p_i[i]
f_b = grid$f_b[i]; f_i = grid$f_i[i]
for (sim in 1:n.sim) {
  print(paste("sim", sim))
  
  if (dgp %in% c("fcomplex", "fLNL")){
    times = 0.2
    data = generate_tutorial_survival_data(n = n, p = p, p_b = p_b, p_i = p_i, f_b = f_b, f_i = f_i, 
                                           dgp = dgp, n.mc = 10, times = times)
    data.test = generate_tutorial_survival_data(n = n.test, p = p, p_b = p_b, p_i = p_i, f_b = f_b, f_i = f_i,
                                                dgp = dgp, n.mc = n.mc, times = times)
  }else if (dgp %in% c("RLsurv1", "RLsurv2", "RLsurv3")){
    if(dgp == "RLsurv1"){
      times = 3.5
    }else if(dgp == "RLsurv2"){
      times = 1
    }else if(dgp == "RLsurv3"){
      times = 1
    }
    data = generate_R_learner_survival_data(n = n, p = p, dgp = dgp, n.mc = 10, times = times)
    data.test = generate_R_learner_survival_data(n = n.test, p = p, dgp = dgp, n.mc = n.mc, times = times)
  }else{
    data = generate_sprint_survival_data(n = n, p = p, dgp = dgp, n.mc = 10)
    data$X = data$X[,-17]         # remove "Diabetes" variable as it is 0 for everyone
    data.test = generate_sprint_survival_data(n = n.test, p = p, dgp = dgp, n.mc = 1)
    data.test$X = data.test$X[,-17]
  }
  
  data$Y = pmax(rep(0.001, length(data$Y)), data$Y)
  true.catesp = data.test$catesp
  true.catesp.sign = data.test$catesp.sign
  
  predictions = matrix(NA, n.test, length(estimators))
  estimator.output = list()
  for (j in 1:length(estimators)) {
    estimator.name = names(estimators)[j]
    print(estimator.name)
    if (dgp %in% c("fcomplex", "fLNL", "CSF2a1", "CSF2b")){
      if(dgp == "fcomplex" | dgp == "fLNL"){
        times = 0.2
        predictions[,j] = as.numeric(unlist(estimators[[estimator.name]](data, data.test, ps = mean(as.numeric(data$W)),
                                                                         times = times, meta_learner = TRUE)))
      }else if(dgp == "type2"){
        times = 0.15
        predictions[,j] = as.numeric(unlist(estimators[[estimator.name]](data, data.test, times = times, meta_learner = TRUE)))
      }else if(dgp == "CSF2b"){
        times = 1
        predictions[,j] = as.numeric(unlist(estimators[[estimator.name]](data, data.test, times = times, meta_learner = TRUE)))
      }
      correct.classification = sign(predictions[,j]) == true.catesp.sign
    }else if (dgp %in% c("RLsurv1", "RLsurv2", "RLsurv3")){
      if(dgp == "RLsurv2"){
        times = 1
        predictions[,j] = as.numeric(unlist(estimators[[estimator.name]](data, data.test, ps = mean(as.numeric(data$W)),
                                                                         times = times, meta_learner = TRUE)))
      }else if(dgp == "RLsurv1"){
        times = 3.5
        predictions[,j] = as.numeric(unlist(estimators[[estimator.name]](data, data.test, times = times, meta_learner = TRUE)))
      }else if(dgp == "RLsurv3"){
        times = 1
        predictions[,j] = as.numeric(unlist(estimators[[estimator.name]](data, data.test, times = times, meta_learner = TRUE)))
      }
      correct.classification = sign(predictions[,j]) == true.catesp.sign
    }else{
      if(dgp %in% c("fixHR20", "hteHR20", "fixHR50", "hteHR50")){
        alpha = 0.05
      }else{
        alpha = 0.01                  # use a lower alpha in (causal_)survival_forest for low event rates (5-10%)
      }
      predictions[,j] = as.numeric(unlist(estimators[[estimator.name]](data, data.test, ps = mean(as.numeric(data$W)),
                                                                       times = 365.25*3, alpha = alpha, meta_learner = TRUE)))
      correct.classification = (2*as.numeric(predictions[,j] > 0.007) - 1) == true.catesp.sign
    }
    
    # calibration slope
    calib_fit <- lm(predictions[,j] ~ true.catesp)
    
    dfj = data.frame(
      estimator.name = estimator.name,
      mse = mean((predictions[,j] - true.catesp)^2),
      bias = mean(abs(predictions[,j] - true.catesp)),
      rcorr = cor(predictions[,j], true.catesp), 
      calib_coef = calib_fit$coefficients[2], 
      concordance = ccdf(predictions[,j], true.catesp), 
      classif.rate = mean(correct.classification, na.rm = TRUE) # NA: to ignore X1 < 0.3 in DGP 4.
    )
    dfj$rcorr[is.na(dfj$rcorr)==TRUE] <- 0  # assign correlation to 0 when CATE = ATE
    estimator.output[[j]] = dfj
  }
  
  # plot the first replicate result
  if (sim==1){
    # if (n==5000 & p==25){
    #   # Density plot of true CATE
    #   png(paste0("Density_true_cate_",dgp, ".png"), width = 5, height = 5, units = 'in', res = 300)
    #   plot(density(true.catesp), xlab = "True CATE", main = dgp)
    #   text(quantile(true.catesp, prob=c(0.8)), 5, paste0("sd = ",round(sd(true.catesp),3)))
    #   dev.off()
    # }
    # Scatter plot of pred and true CATEs
    newnames <- str_replace_all(names(estimators), "estimate_", "") 
    newnames <- str_replace_all(newnames, "ipcw_", "")
    names(estimators) <- names(estimators)
    png(paste0("grid", i, "_", dgp, "_p", p_b, p_i, "_f", f_b, f_i, ".png"), 
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
out.df = do.call(rbind, out)
write.csv(out.df, gzfile(paste0("./ML_grid_",i,"_fcomplex_p1.csv.gz")), row.names = FALSE)
