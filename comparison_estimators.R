# *** Comparison methods ***
# Using modified causal_survival_forest simulation code in "grf" package with the following modifications:
# 1. dgps.R   -- "generate_causal_survival_data" function: the ground truths "cate" and "cate.sign"
#                 are defined as difference in survival probabilities for all five settings
# 2. comparison_estimators.R -- added more estimators, e.g., R-, X-, S-, T-learner, PTO, causal BART
base_surv <- function(fit, Y, D, x, lambda){
  data <- data.frame(t_event=Y, event=D, x)
  tab <- data.frame(table(data[data$event == 1, "t_event"]))
  y <- as.numeric(as.character(sort(unique(tab[,1]))))
  d <- tab[,2]  # number of events at each unique time

  betaHat <- as.vector((fit$glmnet.fit$beta)[,fit$lambda==lambda])
  h0 <- rep(NA, length(y))
  for(l in 1:length(y)){
    h0[l] <- d[l] / sum(exp(x[data$t_event >= y[l], rownames(fit$glmnet.fit$beta)] %*% betaHat))
  }

  S0 <- exp(-cumsum(h0))
  outcome <- data.frame(time=y,survival=S0)
  outcome
}
pred_surv <- function(fit, S0, x, times, lambda){
  link <- predict(fit$glmnet.fit,x,type = "link")[,fit$lambda==lambda]
  colnames(link) <- NULL

  if(length(times)>1){
    S0_t <- rep(NA, length(times))
    for (i in 1:length(times)){
      S0_t[i] <- S0$survival[S0$time>=times[i]][1]
    }
  }else{
    S0_t <- S0$survival[S0$time>=times][1]
  }

  surv <- S0_t^exp(link)
  surv
}
pred_surv_preval <- function(fit, S0, times, lambda){
  link <- fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == lambda]
  colnames(link) <- NULL

  if(length(times)>1){
    S0_t <- rep(NA, length(times))
    for (i in 1:length(times)){
      S0_t[i] <- S0$survival[S0$time>=times[i]][1]
    }
  }else{
    S0_t <- S0$survival[S0$time>=times][1]
  }

  surv <- S0_t^exp(link)
  surv
}

# S-learner
estimate_coxph_sl <- function(data, data.test, times){

  scoxph_fit <- scoxph(x = data$X,
                       w = data$W,
                       y = data$Y,
                       D = data$D,
                       times = times)

  pred_S_coxph <- predict(scoxph_fit, newx = data.test$X, times = times)
  pred_S_coxph
}

estimate_lasso_sl <- function(data, data.test, times){

  slasso_fit <- slasso_surv(x = data$X,
                            w = data$W,
                            y = data$Y,
                            D = data$D,
                            times = times)

  pred_S_lasso <- predict(slasso_fit, newx = data.test$X, times = times)
  pred_S_lasso
}

estimate_grf_sl <- function(data, data.test, times, alpha = 0.05){
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  index <- findInterval(times, Y.grid)
  grffit <- survival_forest(cbind(data$W, data$X),
                            traindat$Y,
                            traindat$D,
                            alpha = alpha,
                            prediction.type = "Nelson-Aalen",
                            failure.times = Y.grid)
  surf1 <- predict(grffit, cbind(rep(1, length(data.test$Y)), data.test$X))$predictions[, index]
  surf0 <- predict(grffit, cbind(rep(0, length(data.test$Y)), data.test$X))$predictions[, index]
  pred_S_grf <- surf1 - surf0
  pred_S_grf
}

# T-learner
estimate_coxph_tl <- function(data, data.test, times){

  traindat <- data.frame(Y = data$Y, D = data$D, W = data$W, data$X)
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]

  # Model for W = 1
  coxph_fit1 <- coxph(Surv(Y, D) ~., data = traindat1)
  bh_dat <- basehaz(coxph_fit1, centered = FALSE)
  bh <- bh_dat[which.min(abs(bh_dat$time - times)),]$hazard
  est_r1 <- predict(coxph_fit1, newdata = data.frame(data.test$X), type="risk")
  surf1 <- exp(-bh)^est_r1

  # Model for W = 0
  coxph_fit0 <- coxph(Surv(Y, D) ~., data = traindat0)
  bh_dat <- basehaz(coxph_fit0, centered = FALSE)
  bh <- bh_dat[which.min(abs(bh_dat$time - times)),]$hazard
  est_r0 <- predict(coxph_fit0, newdata = data.frame(data.test$X), type="risk")
  surf0 <- exp(-bh)^est_r0

  pred_T_coxph <- surf1 - surf0
  pred_T_coxph
}

estimate_lasso_tl <- function(data, data.test, times){

  traindat <- data.frame(Y = data$Y, D = data$D, W = data$W, data$X)
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]

  # Model for W = 1
  foldid <- sample(rep(seq(10), length = length(traindat1$Y)))
  lasso_fit1 <- glmnet::cv.glmnet(as.matrix(traindat1[,3:dim(traindat1)[2]]),
                          Surv(traindat1$Y, traindat1$D),
                          family = "cox",
                          alpha = 1,
                          foldid = foldid)

  S0 <- base_surv(fit = lasso_fit1,
                    Y = traindat1$Y,
                    D = traindat1$D,
                    x = as.matrix(traindat1[,3:dim(traindat1)[2]]),
                    lambda = lasso_fit1$lambda.min)

  surf1 <- pred_surv(fit = lasso_fit1,
                       S0 = S0,
                       x = data.test$X,
                       times = times,
                       lambda = lasso_fit1$lambda.min)


  # Model for W = 0
  foldid <- sample(rep(seq(10), length = length(traindat0$Y)))
  lasso_fit0 <- glmnet::cv.glmnet(as.matrix(traindat0[,3:dim(traindat0)[2]]),
                          Surv(traindat0$Y, traindat0$D),
                          family = "cox",
                          alpha = 1,
                          foldid = foldid)

  S0 <- base_surv(fit = lasso_fit0,
                    Y = traindat0$Y,
                    D = traindat0$D,
                    x = as.matrix(traindat0[,3:dim(traindat0)[2]]),
                    lambda = lasso_fit0$lambda.min)

  surf0 <- pred_surv(fit = lasso_fit0,
                       S0 = S0,
                       x = data.test$X,
                       times = times,
                       lambda = lasso_fit0$lambda.min)

  pred_T_lasso <- surf1 - surf0
  pred_T_lasso
}

estimate_grf_tl <- function(data, data.test, times, alpha = 0.05){
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  index <- findInterval(times, Y.grid)

  # Model for W = 1
  grffit1 <- grf::survival_forest(data$X[dataW==1],
                                  data$Y[dataW==1],
                                  data$D[dataW==1],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen",
                                  failure.times = Y.grid)
  surf1 <- predict(grffit1, data.test$X)$predictions[, index]

  # Model for W = 0
  grffit0 <- grf::survival_forest(data$X[dataW==0],
                                  data$Y[dataW==0],
                                  data$D[dataW==0],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen",
                                  failure.times = Y.grid)
  surf0 <- predict(grffit0, data.test$X)$predictions[, index]

  pred_T_grf <- surf1 - surf0
  pred_T_grf
}


# X-learner lasso
estimate_ipcw_lasso_xl <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){

  # fit model on W==1 (cross-fitted using 'preval' in glmnet)
  foldid1 <- sample(rep(seq(10), length = length(data$Y[dataW==1])))
  lasso_fit1 <- glmnet::cv.glmnet(data$X[dataW==1,],
                                  Surv(data$Y[dataW==1], data$D[dataW==1]),
                                  family = "cox",
                                  alpha = 1,
                                  keep = TRUE,
                                  foldid = foldid1)
  lambda_1_min <- lasso_fit1$lambda[which.min(lasso_fit1$cvm[!is.na(colSums(lasso_fit1$fit.preval))])]
  S0 <- base_surv(lasso_fit1, data$Y[dataW==1], data$D[dataW==1], data$X[dataW==1,], lambda = lambda_1_min)
  surf1 <- rep(NA, length(data$W))
  surf1[data$W==1] <- pred_surv_preval(lasso_fit1, S0, times = times, lambda = lambda_1_min)
  surf1[data$W==0] <- pred_surv(fit = lasso_fit1,
                                S0 = S0,
                                x = data$X[dataW==0,],
                                times = times,
                                lambda = lasso_fit1$lambda.min)

  # fit model on W==0 (cross-fitted using 'preval' in glmnet)
  foldid0 <- sample(rep(seq(10), length = length(data$Y[dataW==0])))
  lasso_fit0 <- glmnet::cv.glmnet(data$X[dataW==0,],
                                  Surv(data$Y[dataW==0], data$D[dataW==0]),
                                  family = "cox",
                                  alpha = 1,
                                  keep = TRUE,
                                  foldid = foldid0)
  lambda_0_min <- lasso_fit0$lambda[which.min(lasso_fit0$cvm[!is.na(colSums(lasso_fit0$fit.preval))])]
  S0 <- base_surv(lasso_fit0, data$Y[dataW==0], data$D[dataW==0], data$X[dataW==0,], lambda = lambda_0_min)
  surf0 <- rep(NA, length(data$W))
  surf0[data$W==0] <- pred_surv_preval(lasso_fit0, S0, times = times, lambda = lambda_0_min)
  surf0[data$W==1] <- pred_surv(fit = lasso_fit0,
                                    S0 = S0,
                                    x = data$X[dataW==1,],
                                    times = times,
                                    lambda = lasso_fit0$lambda.min)

  Tlasso1 <- 1 - surf1
  Tlasso0 <- 1 - surf0

  # IPCW weights (cross-fitted using 'preval' in glmnet)
  if(cen_fit == "KM"){
    shuffle <- sample(length(data$Y))
    kmdat <- data.frame(Y = data$Y[shuffle], D = data$D[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
    c_hat <- rep(NA, nrow(kmdat))
    for(z in 1:10){
      testIndexes <- which(folds==z, arr.ind=TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      c_fit <- survfit(Surv(trainData$Y, 1 - trainData$D) ~ 1)
      cent <- testData$Y
      cent[testData$D==0] <- times
      c_hat[testIndexes] <- summary(c_fit, times = cent)$surv
    }
    shudat <- data.frame(shuffle, c_hat)
    c_hat <- shudat[order(shuffle), ]$c_hat
  }else if (cen_fit == "survival.forest"){
    Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
    c_fit <- grf::survival_forest(cbind(data$W, data$X),
                                  data$Y,
                                  1 - data$D,
                                  failure.times = Y.grid,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
    C.hat <- predict(c_fit)$predictions
    cent <- data$Y; cent[data$D==0] <- times
    cen.times.index <- findInterval(cent, c_fit$failure.times)
    c_hat <- C.hat[cbind(1:length(data$Y), cen.times.index)]
  }

  # Propensity score (cross-fitted using 'preval' in glmnet)
  if (is.null(ps)){
    stop("propensity score needs to be supplied")
  }else{
    ps.train <- rep(ps, length(data$W))
    ps.test <- rep(ps, length(data.test$W))
  }

  weight <- (1 / c_hat) * (1 / ps.train)

  # X-learner
  tempdat <- data.frame(Y = data$Y, D = data$D, W = data$W, weight, X = data$X, Tlasso0, Tlasso1)
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0
  binary_data <- binary_data[complete.cases(binary_data), ]
  b_data <- list(Y = binary_data$Y, D = binary_data$D, W = binary_data$W, X = binary_data$X,
                 wt = binary_data$weight, mu0 = binary_data$Tlasso0, mu1 = binary_data$Tlasso1)

  foldid <- sample(rep(seq(10), length = length(b_data$Y[b_data$W==1])))
  XLfit1 <- glmnet::cv.glmnet(b_data$X[b_data$W==1,],
                              b_data$D[b_data$W==1] - b_data$mu0[b_data$W==1],
                              weights = b_data$wt[b_data$W==1],
                              foldid = foldid,
                              alpha = 1)
  XLtau1 <- as.vector(-predict(XLfit1, data.test$X, s = "lambda.min"))

  foldid <- sample(rep(seq(10), length = length(b_data$Y[b_data$W==0])))
  XLfit0 <- glmnet::cv.glmnet(b_data$X[b_data$W==0,],
                              b_data$mu1[b_data$W==0] - b_data$D[b_data$W==0],
                              weights = b_data$wt[b_data$W==0],
                              foldid = foldid,
                              alpha = 1)
  XLtau0 <- as.vector(-predict(XLfit0, data.test$X, s = "lambda.min"))

  # weighted CATE
  pred_X_lasso <- XLtau1 * (1 - ps.test) + XLtau0 * ps.test
  as.vector(pred_X_lasso)
}

estimate_ipcw_grf_xl <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){

  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  times.index <- findInterval(times, Y.grid)

  # fit model on W==1
  grffit1 <- grf::survival_forest(data$X[data$W==1,],
                                  data$Y[data$W==1],
                                  data$D[data$W==1],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen",
                                  failure.times = Y.grid)
  surf1 <- rep(NA, length(data$W))
  surf1[data$W==1] <- predict(grffit1)$predictions[, times.index]
  surf1[data$W==0] <- predict(grffit1, data$X[data$W==0,])$predictions[, times.index]

  # fit model on W==0
  grffit0 <- grf::survival_forest(data$X[data$W==0,],
                                  data$Y[data$W==0],
                                  data$D[data$W==0],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen",
                                  failure.times = Y.grid)
  surf0 <- rep(NA, length(data$W))
  surf0[traindat$W==0] <- predict(grffit0)$predictions[, times.index]
  surf0[traindat$W==1] <- predict(grffit0, data$X[data$W==1,])$predictions[, times.index]

  Tgrf1 <- 1-surf1
  Tgrf0 <- 1-surf0

  # IPCW weights
  if(cen_fit == "KM"){
    shuffle <- sample(nrow(traindat))
    kmdat <- data.frame(Y = data$Y[shuffle], D = data$D[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks=10, labels=FALSE)
    C.Y.hat <- rep(NA, nrow(kmdat))
    for(z in 1:10){
      testIndexes <- which(folds==z, arr.ind=TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      km_fit <- survfit(Surv(trainData$Y, 1 - trainData$D) ~ 1)
      cent <- testData$Y
      cent[testData$D==0] <- times
      C.Y.hat[testIndexes] <- summary(km_fit, times = cent)$surv
    }
    shudat <- data.frame(shuffle, C.Y.hat)
    C.Y.hat <- shudat[order(shuffle), ]$C.Y.hat
  }else if (cen_fit == "survival.forest"){
    Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
    c_fit <- grf::survival_forest(cbind(data$W, data$X),
                                  data$Y,
                                  1 - data$D,
                                  failure.times = Y.grid,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
    C.hat <- predict(c_fit)$predictions
    cent <- data$Y; cent[data$D==0] <- times
    cen.times.index <- findInterval(cent, c_fit$failure.times)
    C.Y.hat <- C.hat[cbind(1:length(data$Y), cen.times.index)]
  }
  ipcw <- 1 / C.Y.hat

  # Propensity score
  if (is.null(ps)){
    stop("propensity score needs to be supplied")
  }else{
    ps.train <- rep(ps, length(data$W))
    ps_test <- rep(ps, length(data.test$W))
  }

  weight <- ipcw/ps.train  # censoring weight * treatment weight

  # X-learner
  tempdat <- data.frame(Y = data$Y, D = data$D, W = data$W, weight, X = data$X, Tgrf0, Tgrf1)
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0
  binary_data <- binary_data[complete.cases(binary_data), ]
  b_data <- list(Y = binary_data$Y, D = binary_data$D, W = binary_data$W, X = binary_data$X,
                 wt = binary_data$weight, mu0 = binary_data$Tgrf0, mu1 = binary_data$Tgrf1)

  XLfit1 <- grf::regression_forest(b_data$X[b_data$W==1, ],
                                   b_data$D[b_data$W==1] - b_data$mu0[b_data$W==1],
                                   sample.weights = b_data$wt[b_data$W==1])
  XLtau1 <- -predict(XLfit1, data.frame(data.test$X))

  XLfit0 <- grf::regression_forest(b_data$X[b_data$W==0, ],
                                   b_data$mu1[b_data$W==0] - b_data$D[b_data$W==0],
                                   sample.weights = b_data$wt[b_data$W==0])
  XLtau0 <- -predict(XLfit0, data.frame(data.test$X))

  # weighted CATE
  pred_X_grf <- XLtau1 * (1 - ps_test) + XLtau0 * ps_test
  as.vector(pred_X_grf)
}

estimate_ipcw_las_grf_xl <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){
Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
times.index <- findInterval(times, Y.grid)

# fit model on W==1
grffit1 <- grf::survival_forest(data$X[data$W==1,],
                                data$Y[data$W==1],
                                data$D[data$W==1],
                                alpha = alpha,
                                prediction.type = "Nelson-Aalen",
                                failure.times = Y.grid)
surf1 <- rep(NA, length(data$W))
surf1[data$W==1] <- predict(grffit1)$predictions[, times.index]
surf1[data$W==0] <- predict(grffit1, data$X[data$W==0,])$predictions[, times.index]

# fit model on W==0
grffit0 <- grf::survival_forest(data$X[data$W==0,],
                                data$Y[data$W==0],
                                data$D[data$W==0],
                                alpha = alpha,
                                prediction.type = "Nelson-Aalen",
                                failure.times = Y.grid)
surf0 <- rep(NA, length(data$W))
surf0[traindat$W==0] <- predict(grffit0)$predictions[, times.index]
surf0[traindat$W==1] <- predict(grffit0, data$X[data$W==1,])$predictions[, times.index]

Tgrf1 <- 1-surf1
Tgrf0 <- 1-surf0

# IPCW weights
if(cen_fit == "KM"){
  shuffle <- sample(nrow(traindat))
  kmdat <- data.frame(Y = data$Y[shuffle], D = data$D[shuffle])
  folds <- cut(seq(1, nrow(kmdat)), breaks=10, labels=FALSE)
  C.Y.hat <- rep(NA, nrow(kmdat))
  for(z in 1:10){
    testIndexes <- which(folds==z, arr.ind=TRUE)
    testData <- kmdat[testIndexes, ]
    trainData <- kmdat[-testIndexes, ]
    km_fit <- survfit(Surv(trainData$Y, 1 - trainData$D) ~ 1)
    cent <- testData$Y
    cent[testData$D==0] <- times
    C.Y.hat[testIndexes] <- summary(km_fit, times = cent)$surv
  }
  shudat <- data.frame(shuffle, C.Y.hat)
  C.Y.hat <- shudat[order(shuffle), ]$C.Y.hat
}else if (cen_fit == "survival.forest"){
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  c_fit <- grf::survival_forest(cbind(data$W, data$X),
                                data$Y,
                                1 - data$D,
                                failure.times = Y.grid,
                                alpha = alpha,
                                prediction.type = "Nelson-Aalen")
  C.hat <- predict(c_fit)$predictions
  cent <- data$Y; cent[data$D==0] <- times
  cen.times.index <- findInterval(cent, c_fit$failure.times)
  C.Y.hat <- C.hat[cbind(1:length(data$Y), cen.times.index)]
}
ipcw <- 1 / C.Y.hat

# Propensity score
if (is.null(ps)){
  stop("propensity score needs to be supplied")
}else{
  ps.train <- rep(ps, length(data$W))
  ps_test <- rep(ps, length(data.test$W))
}

weight <- ipcw/ps.train  # censoring weight * treatment weight

# X-learner
tempdat <- data.frame(Y = data$Y, D = data$D, W = data$W, weight, X = data$X, Tgrf0, Tgrf1)
binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]
binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0
binary_data <- binary_data[complete.cases(binary_data), ]
b_data <- list(Y = binary_data$Y, D = binary_data$D, W = binary_data$W, X = binary_data$X,
               wt = binary_data$weight, mu0 = binary_data$Tgrf0, mu1 = binary_data$Tgrf1)

foldid <- sample(rep(seq(10), length = length(b_data$Y[b_data$W==1])))
XLfit1 <- glmnet::cv.glmnet(b_data$X[b_data$W==1, ],
                            b_data$D[b_data$W==1] - b_data$mu0[b_data$W==1],
                            weights = b_data$wt[b_data$W==1],
                            foldid = foldid,
                            alpha = 1)
XLtau1 <- as.vector(-predict(XLfit1, data.test$X, s = "lambda.min"))

foldid <- sample(rep(seq(10), length = length(b_data$Y[b_data$W==0])))
XLfit0 <- glmnet::cv.glmnet(b_data$X[b_data$W==0, ],
                            b_data$mu1[b_data$W==0] - b_data$D[b_data$W==0],
                            weights = b_data$wt[b_data$W==0],
                            foldid = foldid,
                            alpha = 1)
XLtau0 <- as.vector(-predict(XLfit0, data.test$X, s = "lambda.min"))

# weighted CATE
pred_X_las_grf <- XLtau1 * (1 - ps_test) + XLtau0 * ps_test
as.vector(pred_X_las_grf)
}


# "IPCW" R-learner
estimate_ipcw_lasso_rl <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){
  rlasso_fit <- rlasso(x = data$X,
                       w = data$W,
                       y = data$Y,
                       D = data$D,
                       p_hat = ps,
                       alpha = alpha,
                       times = times,
                       cen_fit = cen_fit)
  rlasso_est <- predict(object = rlasso_fit, newx = data.test$X)
  as.vector(rlasso_est)
}

estimate_ipcw_grf_rl <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  rgrf_fit <- rgrf(x = data$X,
                   w = data$W,
                   y = data$Y,
                   D = data$D,
                   p_hat = ps,
                   alpha = alpha,
                   failure.times = Y.grid,
                   times = times,
                   cen_fit = cen_fit)
  rgrf_est <- predict(object = rgrf_fit, data.test$X)
  rgrf_est
}

estimate_ipcw_las_grf_rl <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  rlasgrf_fit <- rlasgrf(x = data$X,
                         w = data$W,
                         y = data$Y,
                         D = data$D,
                         p_hat = ps,
                         alpha = alpha,
                         failure.times = Y.grid,
                         times = times,
                         cen_fit = cen_fit)
  rlasgrf_est <- predict(object = rlasgrf_fit, newx = data.test$X)
  rlasgrf_est
}


# "IPCW" F-learner
estimate_ipcw_lasso_fl <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){

  # IPCW weights
  if(cen_fit == "KM"){
    shuffle <- sample(length(data$Y))
    kmdat <- data.frame(Y = data$Y[shuffle], D = data$D[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
    c_hat <- rep(NA, nrow(kmdat))
    for(z in 1:10){
      testIndexes <- which(folds==z, arr.ind=TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      c_fit <- survfit(Surv(trainData$Y, 1 - trainData$D) ~ 1)
      cent <- testData$Y
      cent[testData$D==0] <- times
      c_hat[testIndexes] <- summary(c_fit, times = cent)$surv
    }
    shudat <- data.frame(shuffle, c_hat)
    c_hat <- shudat[order(shuffle), ]$c_hat
  }else if (cen_fit == "survival.forest"){
    Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
    c_fit <- grf::survival_forest(cbind(data$W, data$X),
                                  data$Y,
                                  1 - data$D,
                                  failure.times = Y.grid,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
    C.hat <- predict(c_fit)$predictions
    cent <- data$Y; cent[data$D==0] <- times
    cen.times.index <- findInterval(cent, c_fit$failure.times)
    c_hat <- C.hat[cbind(1:length(data$Y), cen.times.index)]
  }
  ipcw <- 1 / c_hat

  # Propensity score
  if (is.null(ps)){
    stop("propensity score needs to be supplied")
  }else{
    ps_score <- rep(ps, length(data$Y))
  }

  # Subset of uncensored subjects
  tempdat <- data.frame(Y = data$Y, D = data$D, W = data$W, ps_score, ipcw, X = data$X)
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0
  binary_data <- binary_data[complete.cases(binary_data), ]
  b_data <- list(Y = binary_data$Y, D = binary_data$D, W = binary_data$W, X = binary_data$X,
                 wt = binary_data$ipcw, ps = binary_data$ps_score)

  flasso_fit <- Flasso(x = b_data$X,
                       tx = b_data$W,
                       y = b_data$D,
                       pscore = b_data$ps,
                       weight = b_data$wt)
  pred_flasso <- as.vector(-predict(flasso_fit, data.test$X))
  pred_flasso
}

estimate_ipcw_grf_fl <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){

  # IPCW weights
  if(cen_fit == "KM"){
    shuffle <- sample(nrow(traindat))
    kmdat <- data.frame(Y = data$Y[shuffle], D = data$D[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks=10, labels=FALSE)
    C.Y.hat <- rep(NA, nrow(kmdat))
    for(z in 1:10){
      testIndexes <- which(folds==z, arr.ind=TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      km_fit <- survfit(Surv(trainData$Y, 1 - trainData$D) ~ 1)
      cent <- testData$Y
      cent[testData$D==0] <- times
      C.Y.hat[testIndexes] <- summary(km_fit, times = cent)$surv
    }
    shudat <- data.frame(shuffle, C.Y.hat)
    C.Y.hat <- shudat[order(shuffle), ]$C.Y.hat
  }else if (cen_fit == "survival.forest"){
    Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
    c_fit <- grf::survival_forest(cbind(data$W, data$X),
                                  data$Y,
                                  1 - data$D,
                                  failure.times = Y.grid,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
    C.hat <- predict(c_fit)$predictions
    cent <- data$Y; cent[data$D==0] <- times
    cen.times.index <- findInterval(cent, c_fit$failure.times)
    C.Y.hat <- C.hat[cbind(1:length(data$Y), cen.times.index)]
  }
  ipcw <- 1 / C.Y.hat

  # Propensity score
  if (is.null(ps)){
    stop("propensity score needs to be supplied")
  }else{
    pscore <- ps
  }

  # Subset of uncensored subjects
  tempdat <- data.frame(Y = data$Y, D = data$D, W = data$W, pscore, ipcw, X = data$X)
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0
  binary_data <- binary_data[complete.cases(binary_data), ]
  b_data <- list(Y = binary_data$Y, D = binary_data$D, W = binary_data$W, X = binary_data$X,
                 wt = binary_data$ipcw, ps = binary_data$pscore)

  fgrf_fit <- Fgrf(x = b_data$X,
                   tx = b_data$W,
                   y = b_data$D,
                   pscore = b_data$ps,
                   weight = b_data$wt)
  pred_fgrf <- -predict(fgrf_fit, data.frame(data.test$X))
  pred_fgrf
}

# Causal survival forest
estimate_csf_probs <- function(data, data.test, times, alpha = 0.05, ps = NULL, cen_fit = "KM") {
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  fit <- grf::causal_survival_forest(X = data$X,
                                     Y = data$Y,
                                     W = data$W,
                                     D = data$D,
                                     failure.times = Y.grid,
                                     W.hat = ps,
                                     alpha = alpha,
                                     target="survival.probability",
                                     horizon = times)
  p <- predict(fit, data.test$X)
  p$predictions
}
