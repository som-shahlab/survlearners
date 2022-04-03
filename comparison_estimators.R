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
estimate_coxph_sl <- function(data, data.test, times, ps = NULL, alpha = 0.05, cen_fit = "KM", meta_learner = TRUE){
  
  scoxph_fit <- scoxph(x = data$X,
                       w = data$W,
                       y = data$Y, 
                       D = data$D, 
                       times = times)
  
  pred_S_coxph <- predict(scoxph_fit, newx = data.test$X, times = times)
  pred_S_coxph
}

estimate_lasso_sl <- function(data, data.test, times = times, ps = NULL, alpha = 0.05, cen_fit = "KM",
                              lambda_choice = "lambda.min", meta_learner = TRUE){
  
  slasso_fit <- slasso_surv(x = data$X,
                            w = data$W,
                            y = data$Y, 
                            D = data$D,
                            lambda_choice = lambda_choice,
                            times = times)
  
  pred_S_lasso <- predict(slasso_fit, newx = data.test$X, times = times)
  pred_S_lasso
}

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

# T-learner 
estimate_coxph_tl <- function(data, data.test, times, ps = NULL, alpha = 0.05, cen_fit = "KM", meta_learner = TRUE){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
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

estimate_lasso_tl <- function(data, data.test, nfolds = 10, times = times, ps = NULL, alpha = 0.05,
                              cen_fit = "KM", lambda_choice = "lambda.min", meta_learner = TRUE){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  
  # Model for W = 1
  foldid <- sample(rep(seq(nfolds), length = length(traindat1$Y)))
  lasso_fit1 <- cv.glmnet(as.matrix(traindat1[,3:dim(traindat1)[2]]),
                          Surv(traindat1$Y, traindat1$D),
                          family = "cox",
                          alpha = 1,
                          foldid = foldid)
  
  if(lambda_choice == "lambda.min"){
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
  }else if(lambda_choice == "lambda.1se"){
    S0 <- base_surv(fit = lasso_fit1, 
                    Y = traindat1$Y,
                    D = traindat1$D,
                    x = as.matrix(traindat1[,3:dim(traindat1)[2]]),
                    lambda = lasso_fit1$lambda.1se)
    
    surf1 <- pred_surv(fit = lasso_fit1,
                       S0 = S0,
                       x = data.test$X,
                       times = times,
                       lambda = lasso_fit1$lambda.1se)
  }
  
  # Model for W = 0
  foldid <- sample(rep(seq(nfolds), length = length(traindat0$Y)))
  lasso_fit0 <- cv.glmnet(as.matrix(traindat0[,3:dim(traindat0)[2]]),
                          Surv(traindat0$Y, traindat0$D),
                          family = "cox",
                          alpha = 1,
                          foldid = foldid)
  
  if(lambda_choice == "lambda.min"){
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
  }else if(lambda_choice == "lambda.1se"){
    S0 <- base_surv(fit = lasso_fit0, 
                    Y = traindat0$Y,
                    D = traindat0$D,
                    x = as.matrix(traindat0[,3:dim(traindat0)[2]]),
                    lambda = lasso_fit0$lambda.1se)
    
    surf0 <- pred_surv(fit = lasso_fit0,
                       S0 = S0,
                       x = data.test$X,
                       times = times,
                       lambda = lasso_fit0$lambda.1se)
  }
    
  pred_T_lasso <- surf1 - surf0
  pred_T_lasso
}

estimate_grf_tl <- function(data, data.test, times = times, alpha = 0.05, ps = NULL, cen_fit = "KM", meta_learner = TRUE){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  
  Y.grid <- seq(min(traindat$Y), max(traindat$Y), (max(traindat$Y) - min(traindat$Y))/100)
  index <- findInterval(times, Y.grid)
  
  # Model for W = 1
  grffit1 <- survival_forest(as.matrix(traindat1[,3:dim(traindat1)[2]]),
                             traindat1$Y,
                             traindat1$D, 
                             alpha = alpha,
                             prediction.type = "Nelson-Aalen",
                             failure.times = Y.grid)
  surf1 <- predict(grffit1, data.test$X)$predictions[, index]

  # Model for W = 0
  grffit0 <- survival_forest(as.matrix(traindat0[,3:dim(traindat0)[2]]),
                             traindat0$Y,
                             traindat0$D, 
                             alpha = alpha, 
                             prediction.type = "Nelson-Aalen",
                             failure.times = Y.grid)
  surf0 <- predict(grffit0, data.test$X)$predictions[, index]
  
  pred_T_grf <- surf1 - surf0
  pred_T_grf
}


# X-learner lasso
estimate_ipcw_lasso_xl <- function(data, data.test, nfolds = 10, ps = NULL, times = times, alpha = 0.05,
                                   lambda_choice = "lambda.min", cen_fit = "KM", meta_learner = TRUE){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- traindat[, !colnames(traindat) %in% c("W")]   
  
  # fit model on W==1 (cross-fitted using 'preval' in glmnet)
  train1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  foldid1 <- sample(rep(seq(nfolds), length = length(train1$Y)))
  lasso_fit1 <- cv.glmnet(as.matrix(train1[,3:dim(train1)[2]]),
                          Surv(train1$Y, train1$D),
                          family = "cox",
                          alpha = 1,
                          keep = TRUE,
                          foldid = foldid1)

  lambda_1_min <- lasso_fit1$lambda[which.min(lasso_fit1$cvm[!is.na(colSums(lasso_fit1$fit.preval))])]
  S0 <- base_surv(lasso_fit1, train1$Y, train1$D, as.matrix(train1[, 3:dim(train1)[2]]), lambda = lambda_1_min)
  surf1 <- rep(NA, length(traindat$Y))
  surf1[traindat$W==1] <- pred_surv_preval(lasso_fit1, S0, times = times, lambda = lambda_1_min)
  surf1[traindat$W==0] <- pred_surv(fit = lasso_fit1,
                                    S0 = S0,
                                    x = as.matrix(traindat[traindat$W==0, 4:dim(traindat)[2]]),
                                    times = times,
                                    lambda = lasso_fit1$lambda.min)

  # fit model on W==0 (cross-fitted using 'preval' in glmnet)
  train0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  foldid0 <- sample(rep(seq(nfolds), length = length(train0$Y)))
  lasso_fit0 <- cv.glmnet(as.matrix(train0[,3:dim(train0)[2]]),
                          Surv(train0$Y, train0$D),
                          family = "cox",
                          alpha = 1,
                          keep = TRUE,
                          foldid = foldid0)
  
  lambda_0_min <- lasso_fit0$lambda[which.min(lasso_fit0$cvm[!is.na(colSums(lasso_fit0$fit.preval))])]
  S0 <- base_surv(lasso_fit0, train0$Y, train0$D, as.matrix(train0[, 3:dim(train0)[2]]), lambda = lambda_0_min)
  surf0 <- rep(NA, length(traindat$Y))
  surf0[traindat$W==0] <- pred_surv_preval(lasso_fit0, S0, times = times, lambda = lambda_0_min)
  surf0[traindat$W==1] <- pred_surv(fit = lasso_fit0,
                                    S0 = S0,
                                    x = as.matrix(traindat[traindat$W==1, 4:dim(traindat)[2]]),
                                    times = times,
                                    lambda = lasso_fit0$lambda.min)
  
  Tlasso1 <- 1-surf1
  Tlasso0 <- 1-surf0
  
  # IPCW weights (cross-fitted using 'preval' in glmnet)
  if(cen_fit == "KM"){
    shuffle <- sample(nrow(traindat))
    kmdat <- traindat[shuffle,]
    folds <- cut(seq(1, nrow(kmdat)), breaks=10, labels=FALSE)
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
  }else if (cen_fit == "lasso"){
    foldid <- sample(rep(seq(nfolds), length = length(traindat$Y)))
    c_fit <- cv.glmnet(as.matrix(traindat[,3:dim(traindat)[2]]), 
                       Surv(traindat$Y, 1-traindat$D),
                       family = "cox",
                       foldid = foldid,
                       keep = TRUE,
                       alpha = 1)
    c_lambda_min <- c_fit$lambda[which.min(c_fit$cvm[!is.na(colSums(c_fit$fit.preval))])]
    S0 <- base_surv(c_fit, traindat$Y, 1-traindat$D, as.matrix(traindat[,3:dim(traindat)[2]]), lambda = c_lambda_min)
    cent <- rep(times, length(traindat$Y))
    cent[traindat$D==1] <- traindat$Y[traindat$D==1]
    c_hat <- pred_surv_preval(c_fit, S0, times = cent, lambda = c_lambda_min)
  }else if (cen_fit == "survival.forest"){
    Y.grid <- seq(min(traindat$Y), max(traindat$Y), (max(traindat$Y) - min(traindat$Y))/100)
    c_fit <- survival_forest(as.matrix(traindat[,3:dim(traindat)[2]]),
                             traindat$Y,
                             1 - traindat$D,
                             failure.times = Y.grid,
                             alpha = alpha,
                             prediction.type = "Nelson-Aalen")
    C.hat <- predict(c_fit)$predictions
    cent <- traindat$Y
    cent[traindat$D==0] <- times
    cen.times.index <- findInterval(cent, c_fit$failure.times)
    c_hat <- C.hat[cbind(1:length(traindat$Y), cen.times.index)]
  }
  
  # Propensity score (cross-fitted using 'preval' in glmnet)
  if (is.null(ps)){
    foldid <- sample(rep(seq(nfolds), length = length(traindat$Y)))
    w_fit <- cv.glmnet(data$X,
                       traindat$W,
                       family="binomial",
                       type.measure="deviance",
                       foldid = foldid,
                       keep = TRUE,
                       alpha = 1)
    w_lambda_min <- w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
    theta_hat_train <- w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
    ps.train <- 1/(1 + exp(-theta_hat_train)) 
    theta_hat_test <- as.vector(predict(w_fit, newx = data.test$X, s = w_fit$lambda.min))
    ps.test <- 1/(1 + exp(-theta_hat_test))
  }else{
    ps.train <- rep(ps, length(data$W))
    ps.test <- rep(ps, length(data.test$W))
  }
    
  weight <- (1 / c_hat)*(1 / ps.train)
  
  # X-learner
  tempdat <- data.frame(data$Y, data$D, data$W, weight, data$X, Tlasso0, Tlasso1)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]         
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  if (meta_learner){
    foldid <- sample(rep(seq(nfolds), length = nrow(binary_data[binary_data$W==1,])))
    XLfit1 <- cv.glmnet(as.matrix(binary_data[binary_data$W==1, 5:(ncol(binary_data)-2)]), 
                        binary_data$D[binary_data$W==1] - binary_data$Tlasso0[binary_data$W==1],
                        weights = binary_data$weight[binary_data$W==1],
                        foldid = foldid,
                        alpha = 1)
    XLtau1 <- as.vector(-predict(XLfit1, data.test$X, s = lambda_choice))
    
    foldid <- sample(rep(seq(nfolds), length = nrow(binary_data[binary_data$W==0,])))
    XLfit0 <- cv.glmnet(as.matrix(binary_data[binary_data$W==0, 5:(ncol(binary_data)-2)]), 
                        binary_data$Tlasso1[binary_data$W==0] - binary_data$D[binary_data$W==0],
                        weights = binary_data$weight[binary_data$W==0],
                        foldid = foldid,
                        alpha = 1)
    XLtau0 <- as.vector(-predict(XLfit0, data.test$X, s = lambda_choice))
  }else{
    dat1 <- binary_data[binary_data$W==1,]
    datt1 <- data.frame(y1 = dat1$D - dat1$Tlasso0, dat1[, 5:(ncol(dat1)-2)])
    XLfit1 <- glm(y1 ~ ., family = "gaussian", weight = dat1$weight, data = datt1)
    XLtau1 <- -predict(XLfit1, data.frame(data.test$X), meta_learner = meta_learner)
    
    dat0 <- binary_data[binary_data$W==0,]
    datt0 <- data.frame(y0 = dat0$Tlasso1 - dat0$D, dat0[, 5:(ncol(dat0)-2)])
    XLfit0 <- glm(y0 ~ ., family = "gaussian", weight = dat0$weight, data = datt0)
    XLtau0 <- -predict(XLfit0, data.frame(data.test$X), meta_learner = meta_learner)
  }
    
  # weighted CATE
  pred_X_lasso <- XLtau1 * (1 - ps.test) + XLtau0 * ps.test
  as.vector(pred_X_lasso)
}

estimate_ipcw_grf_xl <- function(data, data.test, ps = NULL, times = times, 
                                 alpha = 0.05, cen_fit = "KM", meta_learner = TRUE){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- traindat[, !colnames(traindat) %in% c("W")]   
  
  Y.grid <- seq(min(traindat$Y), max(traindat$Y), (max(traindat$Y) - min(traindat$Y))/100)
  times.index <- findInterval(times, Y.grid)
  
  # fit model on W==1
  train1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  grffit1 <- survival_forest(as.matrix(train1[,3:dim(train1)[2]]),
                             train1$Y,
                             train1$D, 
                             alpha = alpha,
                             prediction.type = "Nelson-Aalen",
                             failure.times = Y.grid)
  surf1 <- rep(NA, length(traindat$Y))
  surf1[traindat$W==1] <- predict(grffit1)$predictions[, times.index]
  surf1[traindat$W==0] <- predict(grffit1, as.matrix(traindat[traindat$W==0, 4:dim(traindat)[2]]))$predictions[, times.index]
  
  # fit model on W==0
  train0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  grffit0 <- survival_forest(as.matrix(train0[,3:dim(train0)[2]]),
                             train0$Y,
                             train0$D,
                             alpha = alpha,
                             prediction.type = "Nelson-Aalen",
                             failure.times = Y.grid)
  surf0 <- rep(NA, length(traindat$Y))
  surf0[traindat$W==0] <- predict(grffit0)$predictions[, times.index]
  surf0[traindat$W==1] <- predict(grffit0, as.matrix(traindat[traindat$W==1, 4:dim(traindat)[2]]))$predictions[, times.index]
  
  Tgbm1 <- 1-surf1
  Tgbm0 <- 1-surf0
  
  # IPCW weights
  if(cen_fit == "KM"){
    shuffle <- sample(nrow(traindat))
    kmdat <- traindat[shuffle,]
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
    sf.censor <- survival_forest(as.matrix(traindat[,3:dim(traindat)[2]]),
                                 traindat$Y,
                                 1 - traindat$D,
                                 failure.times = Y.grid,
                                 alpha = alpha,
                                 prediction.type = "Nelson-Aalen")
    C.hat <- predict(sf.censor)$predictions
    cent <- traindat$Y
    cent[traindat$D==0] <- times
    cen.times.index <- findInterval(cent, sf.censor$failure.times)
    C.Y.hat <- C.hat[cbind(1:length(traindat$Y), cen.times.index)]
  }
  ipcw <- 1 / C.Y.hat
  
  # Propensity score 
  if (is.null(ps)){
    psfit <- regression_forest(data$X, data$W, alpha = alpha) 
    ps.train <- psfit$predictions
    ps_test <- predict(psfit, data.test$X)$predictions
  }else{
    ps.train <- rep(ps, length(data$W))
    ps_test <- rep(ps, length(data.test$W))
  }
  
  weight <- ipcw/ps.train  # censoring weight * treatment weight
    
  # X-learner
  tempdat <- data.frame(data$Y, data$D, data$W, weight, data$X, Tgbm0, Tgbm1)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  if (meta_learner){
    XLfit1 <- regression_forest(as.matrix(binary_data[binary_data$W==1, 5:(ncol(binary_data)-2)]),
                                binary_data$D[binary_data$W==1] - binary_data$Tgbm0[binary_data$W==1],
                                sample.weights = binary_data$weight[binary_data$W==1])
    XLtau1 <- -predict(XLfit1, data.frame(data.test$X))
    
    XLfit0 <- regression_forest(as.matrix(binary_data[binary_data$W==0, 5:(ncol(binary_data)-2)]),
                                binary_data$Tgbm1[binary_data$W==0] - binary_data$D[binary_data$W==0],
                                sample.weights = binary_data$weight[binary_data$W==0])
    XLtau0 <- -predict(XLfit0, data.frame(data.test$X))
  }else{
    dat1 <- binary_data[binary_data$W==1,]
    datt1 <- data.frame(y1 = dat1$D - dat1$Tgbm0, dat1[, 5:(ncol(dat1)-2)])
    XLfit1 <- glm(y1 ~ ., family = "gaussian", weight = dat1$weight, data = datt1)
    XLtau1 <- -predict(XLfit1, data.frame(data.test$X), meta_learner = meta_learner)
    
    dat0 <- binary_data[binary_data$W==0,]
    datt0 <- data.frame(y0 = dat0$Tgbm1 - dat0$D, dat0[, 5:(ncol(dat0)-2)])
    XLfit0 <- glm(y0 ~ ., family = "gaussian", weight = dat0$weight, data = datt0)
    XLtau0 <- -predict(XLfit0, data.frame(data.test$X), meta_learner = meta_learner)
  }
  
  # weighted CATE
  pred_X_grf <- XLtau1 * (1 - ps_test) + XLtau0 * ps_test
  as.vector(pred_X_grf)
}

estimate_ipcw_las_grf_xl <- function(data, data.test, nfolds = 10, ps = NULL, times = times, alpha = 0.05, 
                                     lambda_choice = "lambda.min", cen_fit = "KM", meta_learner = TRUE){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- traindat[, !colnames(traindat) %in% c("W")]   
  
  Y.grid <- seq(min(traindat$Y), max(traindat$Y), (max(traindat$Y) - min(traindat$Y))/100)
  times.index <- findInterval(times, Y.grid)
  
  # fit model on W==1
  train1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  grffit1 <- survival_forest(as.matrix(train1[,3:dim(train1)[2]]),
                             train1$Y,
                             train1$D, 
                             alpha = alpha,
                             prediction.type = "Nelson-Aalen",
                             failure.times = Y.grid)
  surf1 <- rep(NA, length(traindat$Y))
  surf1[traindat$W==1] <- predict(grffit1)$predictions[, times.index]
  surf1[traindat$W==0] <- predict(grffit1, as.matrix(traindat[traindat$W==0, 4:dim(traindat)[2]]))$predictions[, times.index]
  
  # fit model on W==0
  train0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  grffit0 <- survival_forest(as.matrix(train0[,3:dim(train0)[2]]),
                             train0$Y,
                             train0$D,
                             alpha = alpha,
                             prediction.type = "Nelson-Aalen",
                             failure.times = Y.grid)
  surf0 <- rep(NA, length(traindat$Y))
  surf0[traindat$W==0] <- predict(grffit0)$predictions[, times.index]
  surf0[traindat$W==1] <- predict(grffit0, as.matrix(traindat[traindat$W==1, 4:dim(traindat)[2]]))$predictions[, times.index]
  
  Tgbm1 <- 1-surf1
  Tgbm0 <- 1-surf0
  
  # IPCW weights
  if(cen_fit == "KM"){
    shuffle <- sample(nrow(traindat))
    kmdat <- traindat[shuffle,]
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
    sf.censor <- survival_forest(as.matrix(traindat[,3:dim(traindat)[2]]),
                                 traindat$Y,
                                 1 - traindat$D,
                                 failure.times = Y.grid,
                                 alpha = alpha,
                                 prediction.type = "Nelson-Aalen")
    C.hat <- predict(sf.censor)$predictions
    cent <- traindat$Y
    cent[traindat$D==0] <- times
    cen.times.index <- findInterval(cent, sf.censor$failure.times)
    C.Y.hat <- C.hat[cbind(1:length(traindat$Y), cen.times.index)]
  }
  ipcw <- 1 / C.Y.hat
  
  # Propensity score 
  if (is.null(ps)){
    psfit <- regression_forest(data$X, data$W, alpha = alpha) 
    ps.train <- psfit$predictions
    ps_test <- predict(psfit, data.test$X)$predictions
  }else{
    ps.train <- rep(ps, length(data$W))
    ps_test <- rep(ps, length(data.test$W))
  }
  
  weight <- ipcw/ps.train  # censoring weight * treatment weight
  
  # X-learner
  tempdat <- data.frame(data$Y, data$D, data$W, weight, data$X, Tgbm0, Tgbm1)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  if (meta_learner){
    foldid <- sample(rep(seq(nfolds), length = nrow(binary_data[binary_data$W==1,])))
    XLfit1 <- cv.glmnet(as.matrix(binary_data[binary_data$W==1, 5:(ncol(binary_data)-2)]), 
                        binary_data$D[binary_data$W==1] - binary_data$Tgbm0[binary_data$W==1],
                        weights = binary_data$weight[binary_data$W==1],
                        foldid = foldid,
                        alpha = 1)
    XLtau1 <- as.vector(-predict(XLfit1, data.test$X, s = lambda_choice))
    
    foldid <- sample(rep(seq(nfolds), length = length(binary_data[binary_data$W==0,]$Y)))
    XLfit0 <- cv.glmnet(as.matrix(binary_data[binary_data$W==0, 5:(ncol(binary_data)-2)]), 
                        binary_data$Tgbm1[binary_data$W==0] - binary_data$D[binary_data$W==0],
                        weights = binary_data$weight[binary_data$W==0],
                        foldid = foldid,
                        alpha = 1)
    XLtau0 <- as.vector(-predict(XLfit0, data.test$X, s = lambda_choice))
  }else{
    dat1 <- binary_data[binary_data$W==1,]
    datt1 <- data.frame(y1 = dat1$D - dat1$Tgbm0, dat1[, 5:(ncol(dat1)-2)])
    XLfit1 <- glm(y1 ~ ., family = "gaussian", weight = dat1$weight, data = datt1)
    XLtau1 <- -predict(XLfit1, data.frame(data.test$X), meta_learner = meta_learner)
    
    dat0 <- binary_data[binary_data$W==0,]
    datt0 <- data.frame(y0 = dat0$Tgbm1 - dat0$D, dat0[, 5:(ncol(dat0)-2)])
    XLfit0 <- glm(y0 ~ ., family = "gaussian", weight = dat0$weight, data = datt0)
    XLtau0 <- -predict(XLfit0, data.frame(data.test$X), meta_learner = meta_learner)
  }
  
  # weighted CATE
  pred_X_las_grf <- XLtau1 * (1 - ps_test) + XLtau0 * ps_test
  as.vector(pred_X_las_grf)
}


# "IPCW" R-learner 
estimate_ipcw_lasso_rl <- function(data, data.test, ps = NULL, times=times, alpha = 0.05, k_folds = 10, 
                                   lambda_choice = "lambda.min", cen_fit = "KM", meta_learner = TRUE){
  rlasso_fit <- rlasso(x = data$X, 
                       w = data$W, 
                       y = data$Y, 
                       D = data$D, 
                       p_hat = ps,
                       alpha = 1,
                       k_folds = k_folds,
                       times = times, 
                       lambda_choice = lambda_choice, # lambds choice on tau
                       cen_fit = cen_fit,
                       meta_learner = meta_learner)
  rlasso_est <- predict(object = rlasso_fit, newx = data.test$X) 
  as.vector(rlasso_est)
}

estimate_ipcw_grf_rl <- function(data, data.test, ps = NULL, times=times,
                                 alpha = 0.05, cen_fit = "KM", meta_learner = TRUE){
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  rgrf_fit <- rgrf(x = data$X, 
                   w = data$W, 
                   y = data$Y, 
                   D = data$D,
                   p_hat = ps,
                   alpha = alpha, 
                   failure.times = Y.grid,
                   times = times, 
                   cen_fit = cen_fit,
                   meta_learner = meta_learner)
  rgrf_est <- predict(object = rgrf_fit, data.test$X, meta_learner = meta_learner) 
  rgrf_est
}

estimate_ipcw_las_grf_rl <- function(data, data.test, ps = NULL, times=times, k_folds = 10, lambda_choice = "lambda.min",
                                     alpha = 0.05, cen_fit = "KM", meta_learner = TRUE){
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  rlasgrf_fit <- rlasgrf(x = data$X, 
                         w = data$W, 
                         y = data$Y, 
                         D = data$D,
                         p_hat = ps,
                         alpha = alpha, 
                         failure.times = Y.grid,
                         times = times, 
                         k_folds = k_folds,
                         lambda_choice = lambda_choice,
                         cen_fit = cen_fit,
                         meta_learner = meta_learner)
  rlasgrf_est <- predict(object = rlasgrf_fit, newx = data.test$X) 
  rlasgrf_est
}


# "IPCW" F-learner
estimate_ipcw_lasso_fl <- function(data, data.test, nfolds = 10, ps = NULL, times = times, alpha = 0.05,
                                       lambda_choice = "lambda.min", cen_fit = "KM", meta_learner = TRUE){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  foldid <- sample(rep(seq(nfolds), length = length(data$Y)))
  
  # IPCW weights
  if(cen_fit == "KM"){
    shuffle <- sample(nrow(traindat))
    kmdat <- traindat[shuffle,]
    folds <- cut(seq(1, nrow(kmdat)), breaks=10, labels=FALSE)
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
  }else if (cen_fit == "lasso"){
    c_fit <- cv.glmnet(as.matrix(traindat[,3:dim(traindat)[2]]), 
                       Surv(traindat$Y, 1-traindat$D),
                       family = "cox",
                       foldid = foldid,
                       keep = TRUE,
                       alpha = 1)
    
    c_lambda_min <- c_fit$lambda[which.min(c_fit$cvm[!is.na(colSums(c_fit$fit.preval))])]
    S0 <- base_surv(c_fit, traindat$Y, 1-traindat$D, as.matrix(traindat[,3:dim(traindat)[2]]), lambda = c_lambda_min)
    cent <- rep(times, length(traindat$Y))
    cent[traindat$D==1] <- traindat$Y[traindat$D==1]
    c_hat <- pred_surv_preval(c_fit, S0, times = cent, lambda = c_lambda_min)
  }else if (cen_fit == "survival.forest"){
    Y.grid <- seq(min(traindat$Y), max(traindat$Y), (max(traindat$Y) - min(traindat$Y))/100)
    c_fit <- survival_forest(as.matrix(traindat[,3:dim(traindat)[2]]),
                                 traindat$Y,
                                 1 - traindat$D,
                                 failure.times = Y.grid,
                                 prediction.type = "Nelson-Aalen", 
                                 alpha = alpha)
    C.hat <- predict(c_fit)$predictions
    cent <- traindat$Y
    cent[traindat$D==0] <- times
    cen.times.index <- findInterval(cent, Y.grid)
    c_hat <- C.hat[cbind(1:length(traindat$Y), cen.times.index)]
  }
  ipcw <- 1 / c_hat
  
  # Propensity score
  if (is.null(ps)){
    w_fit <- cv.glmnet(data$X, 
                       data$W,
                       family="binomial",
                       type.measure="deviance",
                       foldid = foldid,
                       keep = TRUE,
                       alpha = 1)
    w_lambda_min <- w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
    theta_hat <- w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
    ps_score <- 1/(1 + exp(-theta_hat))
  }else{
    ps_score <- ps
  }
    
  # Subset of uncensored subjects
  tempdat <- data.frame(data$Y, data$D, data$W, ps_score, ipcw, data$X)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0    
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  flasso_fit <- Flasso(x = as.matrix(binary_data[,6:dim(binary_data)[2]]), 
                       tx = binary_data$W, 
                       y = binary_data$D, 
                       pscore = binary_data$ps_score, 
                       weight = binary_data$ipcw, 
                       alpha = 1,
                       nfolds = nfolds, 
                       meta_learner = meta_learner)
  pred_flasso <- as.vector(-predict(flasso_fit, data.test$X, lambda_choice = lambda_choice, meta_learner = meta_learner))
  pred_flasso
}

estimate_ipcw_grf_fl <- function(data, data.test, ps = NULL, times = times, alpha = 0.05,
                                 cen_fit = "KM", meta_learner = TRUE){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")

  # IPCW weights
  if(cen_fit == "KM"){
    shuffle <- sample(nrow(traindat))
    kmdat <- traindat[shuffle,]
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
    Y.grid <- seq(min(traindat$Y), max(traindat$Y), (max(traindat$Y) - min(traindat$Y))/100)
    sf.censor <- survival_forest(as.matrix(traindat[,3:dim(traindat)[2]]),
                                 traindat$Y,
                                 1 - traindat$D,
                                 failure.times = Y.grid,
                                 prediction.type = "Nelson-Aalen", 
                                 alpha = alpha)
    C.hat <- predict(sf.censor)$predictions
    cent <- traindat$Y
    cent[traindat$D==0] <- times
    cen.times.index <- findInterval(cent, Y.grid)
    C.Y.hat <- C.hat[cbind(1:length(traindat$Y), cen.times.index)]
  }
  ipcw <- 1 / C.Y.hat
  
  # Propensity score
  if (is.null(ps)){
    psfit <- regression_forest(as.matrix(traindat[,4:dim(traindat)[2]]),
                               traindat$W, 
                               alpha = alpha) 
    pscore <- psfit$predictions
  }else{
    pscore <- ps
  }
  
  # Subset of uncensored subjects
  tempdat <- data.frame(data$Y, data$D, data$W, pscore, ipcw, data$X)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  fgrf_fit <- Fgrf(x = as.matrix(binary_data[,6:dim(binary_data)[2]]), 
                   tx = binary_data$W, 
                   y = binary_data$D, 
                   pscore = binary_data$pscore, 
                   weight = binary_data$ipcw, 
                   alpha = alpha, 
                   meta_learner = meta_learner)
  pred_fgrf <- -predict(fgrf_fit, data.frame(data.test$X), meta_learner = meta_learner)
  pred_fgrf
}

# Causal survival forest 
estimate_csf_probs <- function(data, data.test, ps = NULL, times = times, alpha = 0.05, cen_fit = "KM", meta_learner = TRUE) {
  Y.grid <- seq(min(data$Y), max(data$Y), (max(data$Y) - min(data$Y))/100)
  fit <- causal_survival_forest(X = data$X,
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

