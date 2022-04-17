rm(list = ls())
library(survlearners)
library(stringr)
source("comparison_estimators.R")

# *** Comparison methods ***
estimators <- list(cate_sl_coxph = cate_sl_coxph,
                   cate_tl_coxph = cate_tl_coxph,
                   cate_csf_probs = cate_csf_probs,
                   cate_xl_grf_lasso = cate_xl_grf_lasso,
                   cate_rl_grf_lasso = cate_rl_grf_lasso,

                   cate_sl_lasso = cate_sl_lasso,
                   cate_tl_lasso = cate_tl_lasso,
                   cate_fl_lasso = cate_fl_lasso,
                   cate_xl_lasso = cate_xl_lasso,
                   cate_rl_lasso = cate_rl_lasso,

                   cate_sl_grf = cate_sl_grf,
                   cate_tl_grf = cate_tl_grf,
                   cate_fl_grf = cate_fl_grf,
                   cate_xl_grf = cate_xl_grf,
                   cate_rl_grf = cate_rl_grf,

                   cate_xl_grf_lasso_sf = cate_xl_grf_lasso_sf,
                   cate_rl_grf_lasso_sf = cate_rl_grf_lasso_sf,
                   cate_fl_lasso_sf = cate_fl_lasso_sf,
                   cate_xl_lasso_sf = cate_xl_lasso_sf,
                   cate_rl_lasso_sf = cate_rl_lasso_sf,
                   cate_fl_grf_sf = cate_fl_grf_sf,
                   cate_xl_grf_sf = cate_xl_grf_sf,
                   cate_rl_grf_sf = cate_rl_grf_sf)

# *** Setup ***
out <- list()
n.sim <- 20
n.mc <- 10000

# Simulations scenarios
grid <- expand.grid(n = 5000,
                    p = 25,
                    n.test = 5000,
                    dgp = c("fcomplex"),
                    p.b = c(1, 25, 25),
                    f.b = c("L", "NL", "NL"),
                    pi = c(0.5),
                    gamma = 1,
                    rho = c(2),
                    cen.scale = c(4),
                    cenM = c("indX"),
                    t0 = (0.2),
                    stringsAsFactors = FALSE)
grid$f.i <- c(rep("L", 6), rep("NL", 3))
grid$p.i <- rep(c(1, 1, 25), 3)
grid <- rbind(grid, grid[2,], grid[5,], grid[8,])
grid[10:12, ]$pi <- c(0.05)                           # unbalanced design
grid <- rbind(grid, grid[1,], grid[1,], grid[1,], grid[4,], grid[4,], grid[4,], grid[4,])
grid[c(13:14, 16:17), ]$cen.scale <- c(8, 7, 8, 7)    # vary censoring rate (under indX): 30% (default), 70% (early censor), 65%
grid[c(13, 16), ]$rho <- rep(1, 2)
grid[c(15, 18), ]$cenM <- rep("dX", 2)                # vary censoring generating model (= dX)
grid[19, ]$cenM <- "dX.ub"                            # unbalanced censoring
grid <- rbind(grid, grid[2,], grid[5,], grid[8,], grid[2,], grid[5,], grid[8,])  # vary heterogeneity: sd(CATE)/sd(mu0sp) 0.17, 0.55, 0.9 (baseline)
grid[20:25, ]$gamma <- c(rep(0.46, 3), rep(0, 3))
grid <- rbind(grid, grid[2,], grid[5,], grid[8,], grid[2,], grid[5,], grid[8,], grid[2,], grid[5,], grid[8,])  # vary event rate
grid[26:34, ]$t0 <- c(rep(0.02,3), rep(0.001,3), rep(0.35,3))
rownames(grid) <- 1:dim(grid)[1]

if(length(args <- commandArgs(T))>0){
  stopifnot(length(args)==1)
  i <- as.integer(args[[1]])
  message("running for grid ", i)
}

n <- grid$n[i]
p <- grid$p[i]
n.test <- grid$n.test[i]
pi <- grid$pi[i]
dgp <- grid$dgp[i]
p.b <- grid$p.b[i]; p.i <- grid$p.i[i]
f.b <- grid$f.b[i]; f.i <- grid$f.i[i]
gamma <- grid$gamma[i]
rho <- grid$rho[i]
cen.scale <- grid$cen.scale[i]
cenM <- grid$cenM[i]
t0 <- grid$t0[i]
an.error.occured <- rep(NA, n.sim)
for (sim in 1:n.sim) {
  #tryCatch( {
  print(paste("sim", sim))
  data <- generate_tutorial_survival_data(n = n, p = p, p.b = p.b, p.i = p.i, f.b = f.b, f.i = f.i, pi = pi,
                                          gamma = gamma, rho = rho, cen.scale = cen.scale, cenM = cenM, dgp = dgp,
                                          n.mc = 10, t0 = t0)
  data.test <- generate_tutorial_survival_data(n = n, p = p, p.b = p.b, p.i = p.i, f.b = f.b, f.i = f.i, pi = pi,
                                               gamma = gamma, rho = rho, cen.scale = cen.scale, cenM = cenM, dgp = dgp,
                                               n.mc = n.mc, t0 = t0)

  data$Y <- pmax(rep(0.001, length(data$Y)), data$Y)
  true.catesp <- data.test$catesp
  true.catesp.sign <- data.test$catesp.sign

  predictions <- matrix(NA, n.test, length(estimators))
  estimator.output <- list()
  for (j in 1:length(estimators)) {
    estimator.name <- names(estimators)[j]
    print(estimator.name)
    if (grepl("sl", estimator.name, fixed = TRUE) == TRUE | grepl("tl", estimator.name, fixed = TRUE) == TRUE) {
      predictions[,j] <- estimators[[estimator.name]](data, data.test, t0 = t0)
    } else {
      predictions[,j] <- estimators[[estimator.name]](data, data.test, t0 = t0, W.hat = pi)
    }

    correct.classification <- sign(predictions[,j]) == true.catesp.sign

    # calibration slope
    calib.fit <- lm(predictions[,j] ~ true.catesp)

    dfj <- data.frame(estimator.name = estimator.name,
                      true.catesp.var = var(true.catesp),
                      mse = mean((predictions[,j] - true.catesp)^2),
                      bias = mean(abs(predictions[,j] - true.catesp)),
                      rcorr = cor(predictions[,j], true.catesp),
                      taucorr = cor(predictions[,j], true.catesp, method = "kendall"),  # use Kendall's tau for concordance
                      calib.coef = calib.fit$coefficients[2],
                      classif.rate = mean(correct.classification, na.rm = TRUE) # NA: to ignore X1 < 0.3 in DGP 4.
    )
    dfj$rcorr[is.na(dfj$rcorr)==TRUE] <- 0  # assign correlation to 0 when CATE = ATE
    estimator.output[[j]] <- dfj
  }

  # Scatter plot of pred and true CATEs
  if (sim==1){
    newnames <- str_replace_all(names(estimators), "cate_", "")
    names(estimators) <- names(estimators)
    png(paste0("grid", i, "cen.fit.png"),
        width = 10, height = 10, units = 'in', res = 300)
    par(mfrow=c(5,5),
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
      axis(side = 1, at=levels, labels = if (z %in% 21:25) levels else FALSE)
      axis(side = 2, at=levels, labels = if (z %in% c(1, 6, 11, 16, 21)) levels else FALSE)
    }
    title(xlab = "Estimated CATE",
          ylab = "True CATE",
          cex.lab = 1.5,
          outer = TRUE, line = 2)
    dev.off()
  }

  df <- do.call(rbind, estimator.output)
  df$n <- n
  df$p <- p
  df$n.test <- n.test
  df$dgp <- dgp
  df$horizon <- t0
  df$sim <- sim

  out <- c(out, list(df))
  #}
  #, error = function(e) {an.error.occured[sim] <<- TRUE})
}
print(sum(an.error.occured, na.rm = TRUE))
out.df <- do.call(rbind, out)
write.csv(out.df, gzfile(paste0("./ML_grid_",i,"_SFfit_p1.csv.gz")), row.names = FALSE)
