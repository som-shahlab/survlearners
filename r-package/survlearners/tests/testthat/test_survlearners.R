library(survlearners)
library(grf)

n = 5000
p = 25
n.test = 5000
p.b = 1
p.i = 1
f.b = "NL"
f.i = "NL"
dgp = "fcomplex"

t0 <- 0.2
data <- survlearners:::generate_tutorial_survival_data(n = n, p = p, p.b = p.b, p.i = p.i, f.b = f.b, f.i = f.i,
                                       dgp = dgp, n.mc = 10, times = t0)

data.test <- survlearners:::generate_tutorial_survival_data(n = n.test, p = p, p.b = p.b, p.i = p.i, f.b = f.b, f.i = f.i,
                                            dgp = dgp, n.mc = 10000, times = t0)


# surv_s_grf: Compare new implementation against Crystal's previous one


estimate_grf_sl <- function(data, data.test, times = times, alpha = 0.05, ps = NULL, cen_fit = "KM", meta_learner = TRUE){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  testdat <- data.frame(data.test$Y, data.test$D, data.test$W, data.test$X)
  colnames(traindat)[1:3] <-colnames(testdat)[1:3] <- c("Y", "D", "W")
  testdat1 <- testdat; testdat1$W <- 1
  testdat0 <- testdat; testdat0$W <- 0

  Y.grid <- seq(min(traindat$Y), max(traindat$Y), (max(traindat$Y) - min(traindat$Y))/100)
  index <- findInterval(times, Y.grid)
  grffit <- survival_forest(as.matrix(traindat[,3:dim(traindat)[2]]),
                            traindat$Y,
                            traindat$D,
                            alpha = alpha,
                            prediction.type = "Nelson-Aalen",
                            failure.times = Y.grid)
  surf1 <- predict(grffit, as.matrix(testdat1[,3:dim(testdat1)[2]]))$predictions[, index]
  surf0 <- predict(grffit, as.matrix(testdat0[,3:dim(testdat0)[2]]))$predictions[, index]
  pred_S_grf <- surf1 - surf0
  pred_S_grf
}

# set.seed(1)
# tau_hat <- surv_s_grf(data, data.test, t0)
set.seed(1)
tau_hat_old <- estimate_grf_sl(data, data.test, times = t0)

# expect_equal(tau_hat, tau_hat_old)
