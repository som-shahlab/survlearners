#' @title F-learner of lasso
#'
#' @description  F-learner, implemented via glmnet (lasso) with 'coxph' distribution
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
#' @param ps The propensity score
#' @param cen.fit The choice of model fitting for censoring
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
#' surv.fl.lasso.fit = surv_fl_lasso(X, W, Y, D, times, ps = 0.5)
#' cate = predict(surv.fl.lasso.fit)
#' cate.test = predict(surv.fl.lasso.fit, X.test)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_fl_lasso <- function(X, W, Y, D, times, alpha = 0.05, ps = NULL, cen.fit = "KM"){

  # IPCW weights
  if(cen.fit == "KM"){
    shuffle <- sample(length(Y))
    kmdat <- data.frame(Y = Y[shuffle], D = D[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
    c.hat <- rep(NA, nrow(kmdat))
    for(z in 1:10){
      testIndexes <- which(folds==z, arr.ind=TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      c.fit <- survival::survfit(survival::Surv(trainData$Y, 1 - trainData$D) ~ 1)
      cent <- testData$Y
      cent[testData$D==0] <- times
      c.hat[testIndexes] <- summary(c.fit, times = cent)$surv
    }
    shudat <- data.frame(shuffle, c.hat)
    c.hat <- shudat[order(shuffle), ]$c.hat
  }else if (cen.fit == "survival.forest"){
    c.fit <- grf::survival_forest(cbind(W, X),
                                  Y,
                                  1 - D,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
    C.hat <- predict(c.fit)$predictions
    cent <- Y; cent[D==0] <- times
    cen.times.index <- findInterval(cent, c.fit$failure.times)
    c.hat <- C.hat[cbind(1:length(Y), cen.times.index)]
  }
  ipcw <- 1 / c.hat

  # Propensity score
  if (is.null(ps)){
    stop("propensity score needs to be supplied")
  }else{
    ps.score <- rep(ps, length(Y))
  }

  # Subset of uncensored subjects
  tempdat <- data.frame(Y = Y, D = D, W = W, ps.score, ipcw, X)
  binary.data <- tempdat[tempdat$D==1|tempdat$Y > times,]
  binary.data$D[binary.data$D==1 & binary.data$Y > times] <- 0
  binary.data <- binary.data[complete.cases(binary.data), ]
  b.data <- list(Y = binary.data$Y, D = binary.data$D, W = binary.data$W,
                 X = as.matrix(binary.data[,6:ncol(binary.data)]),
                 wt = binary.data$ipcw, ps = binary.data$ps.score)

  Z <- b.data$W * b.data$D / b.data$ps - (1 - b.data$W) * b.data$D / (1 - b.data$ps)
  flasso.fit <- glmnet::cv.glmnet(b.data$X, Z, family = "gaussian", weights = b.data$wt, nfolds = 10, alpha = 1)
  flasso.tau <- -predict(flasso.fit, X)

  ret <- list(fit = flasso.fit,
              tau = flasso.tau)
  class(ret) <- 'surv_fl_lasso'
  ret
}

#' predict for surv_fl_lasso
#'
#' get estimated tau(X) using the trained surv_fl_lasso model
#'
#' @param object An surv_fl_lasso object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
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
#' surv.fl.lasso.fit = surv_fl_lasso(X, W, Y, D, times, ps = 0.5)
#' cate = predict(surv.fl.lasso.fit)
#' cate.test = predict(surv.fl.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_fl_lasso = function(object,
                                 newdata = NULL,
                                 ...) {
  if(is.null(newdata)){
    return(object$tau)
  }else{
    return(-predict(object$fit, newdata))
  }
}
