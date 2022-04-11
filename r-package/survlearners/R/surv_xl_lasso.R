#' @title X-learner of lasso
#'
#' @description  X-learner, implemented via glmnet (lasso) with 'coxph' distribution
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
#' @param ps The propensity score
#' @param cen_fit The choice of model fitting for censoring
#' @param newX The test data set (covariates only)
#' @examples
#' \donttest{
#' n = 1000; p = 25
#' times = 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[,1] + (-0.5 - 1 * X[,2]) * W))^2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC/(4^2))^(1/2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv_xl_lasso_fit = surv_xl_lasso(X, W, Y, D, times, ps = 0.5)
#' cate = predict(surv_xl_lasso_fit)
#' cate.test = predict(surv_xl_lasso_fit, X.test)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_xl_lasso <- function(X, W, Y, D, times, alpha = 0.05, ps = NULL, cen_fit = "KM"){

  # fit model on W==1 (cross-fitted using 'preval' in glmnet)
  foldid1 <- sample(rep(seq(10), length = length(Y[W==1])))
  x1 <- as.matrix(data.frame(X[W==1, ]))
  lasso_fit1 <- glmnet::cv.glmnet(x1,
                                  survival::Surv(Y[W==1], D[W==1]),
                                  family = "cox",
                                  alpha = 1,
                                  keep = TRUE,
                                  foldid = foldid1)
  lambda_1_min <- lasso_fit1$lambda[which.min(lasso_fit1$cvm[!is.na(colSums(lasso_fit1$fit.preval))])]
  S0 <- base_surv(lasso_fit1, Y[W==1], D[W==1], x1, lambda = lambda_1_min)
  surf1 <- rep(NA, length(W))
  surf1[W==1] <- pred_surv_preval(lasso_fit1, S0, times = times, lambda = lambda_1_min)
  surf1[W==0] <- pred_surv(fit = lasso_fit1,
                                S0 = S0,
                                X = X[W==0, ],
                                times = times,
                                lambda = lasso_fit1$lambda.min)

  # fit model on W==0 (cross-fitted using 'preval' in glmnet)
  foldid0 <- sample(rep(seq(10), length = length(Y[W==0])))
  x0 <- as.matrix(data.frame(X[W==0, ]))
  lasso_fit0 <- glmnet::cv.glmnet(x0,
                                  survival::Surv(Y[W==0], D[W==0]),
                                  family = "cox",
                                  alpha = 1,
                                  keep = TRUE,
                                  foldid = foldid0)
  lambda_0_min <- lasso_fit0$lambda[which.min(lasso_fit0$cvm[!is.na(colSums(lasso_fit0$fit.preval))])]
  S0 <- base_surv(lasso_fit0, Y[W==0], D[W==0], x0, lambda = lambda_0_min)
  surf0 <- rep(NA, length(W))
  surf0[W==0] <- pred_surv_preval(lasso_fit0, S0, times = times, lambda = lambda_0_min)
  surf0[W==1] <- pred_surv(fit = lasso_fit0,
                                S0 = S0,
                                X = X[W==1, ],
                                times = times,
                                lambda = lasso_fit0$lambda.min)

  Tlasso1 <- 1 - surf1
  Tlasso0 <- 1 - surf0

  # IPCW weights (cross-fitted using 'preval' in glmnet)
  if(cen_fit == "KM"){
    shuffle <- sample(length(Y))
    kmdat <- data.frame(Y = Y[shuffle], D = D[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
    c_hat <- rep(NA, nrow(kmdat))
    for(z in 1:10){
      testIndexes <- which(folds==z, arr.ind=TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      c_fit <- survival::survfit(survival::Surv(trainData$Y, 1 - trainData$D) ~ 1)
      cent <- testData$Y
      cent[testData$D==0] <- times
      c_hat[testIndexes] <- summary(c_fit, times = cent)$surv
    }
    shudat <- data.frame(shuffle, c_hat)
    c_hat <- shudat[order(shuffle), ]$c_hat
  }else if (cen_fit == "survival.forest"){
    c_fit <- grf::survival_forest(cbind(W, X),
                                  Y,
                                  1 - D,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
    C.hat <- predict(c_fit)$predictions
    cent <- Y; cent[D==0] <- times
    cen.times.index <- findInterval(cent, c_fit$failure.times)
    c_hat <- C.hat[cbind(1:length(Y), cen.times.index)]
  }

  # Propensity score (cross-fitted using 'preval' in glmnet)
  if (is.null(ps)){
    stop("propensity score needs to be supplied")
  }else{
    ps.train <- rep(ps, length(W))
  }

  weight <- (1 / c_hat) * (1 / ps.train)

  # X-learner
  tempdat <- data.frame(Y = Y, D = D, W = W, weight, X, Tlasso0, Tlasso1)
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0
  binary_data <- binary_data[complete.cases(binary_data), ]
  b_data <- list(Y = binary_data$Y, D = binary_data$D, W = binary_data$W,
                 X = as.matrix(binary_data[,5:(ncol(binary_data)-2)]),
                 wt = binary_data$weight, mu0 = binary_data$Tlasso0, mu1 = binary_data$Tlasso1)

  foldid <- sample(rep(seq(10), length = length(b_data$Y[b_data$W==1])))
  XLfit1 <- glmnet::cv.glmnet(b_data$X[b_data$W==1, ],
                              b_data$D[b_data$W==1] - b_data$mu0[b_data$W==1],
                              weights = b_data$wt[b_data$W==1],
                              foldid = foldid,
                              alpha = 1)
  XLtau1 <- as.vector(-predict(XLfit1, X, s = "lambda.min"))

  foldid <- sample(rep(seq(10), length = length(b_data$Y[b_data$W==0])))
  XLfit0 <- glmnet::cv.glmnet(b_data$X[b_data$W==0, ],
                              b_data$mu1[b_data$W==0] - b_data$D[b_data$W==0],
                              weights = b_data$wt[b_data$W==0],
                              foldid = foldid,
                              alpha = 1)
  XLtau0 <- as.vector(-predict(XLfit0, X, s = "lambda.min"))

  # weighted CATE
  pred_X_lasso <- as.vector(XLtau1 * (1 - ps.train) + XLtau0 * ps.train)

  ret <- list(fit1 = XLfit1,
              fit0 = XLfit0,
              ps = ps.train,
              tau = pred_X_lasso)
  class(ret) <- 'surv_xl_lasso'
  ret
}

#' predict for surv_xl_lasso
#'
#' get estimated tau(X) using the trained surv_xl_lasso model
#'
#' @param object An surv_xl_lasso object
#' @param newx Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param ps The propensity score
#' @param ... Additional arguments (currently not used)
#'
#' @examples
#' \donttest{
#' n = 1000; p = 25
#' times = 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[,1] + (-0.5 - 1 * X[,2]) * W))^2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC/(4^2))^(1/2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv_xl_lasso_fit = surv_xl_lasso(X, W, Y, D, times, ps = 0.5)
#' cate = predict(surv_xl_lasso_fit)
#' cate.test = predict(surv_xl_lasso_fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_xl_lasso = function(object,
                                 newx = NULL,
                                 ps = NULL,
                                 ...) {
  if(is.null(newx)){
    return(object$tau)
  }else{
    XLtau1 <- as.vector(-predict(object$fit1, newx, s = "lambda.min"))
    XLtau0 <- as.vector(-predict(object$fit0, newx, s = "lambda.min"))
    if(is.null(ps)){
      ps <- rep(object$ps[1], nrow(newx))
    }else{
      if(length(ps)==1){
      ps <- rep(ps, nrow(newx))
      }
    }
    return(as.vector(XLtau1 * (1 - ps) + XLtau0 * ps))
  }
}
