#' @title X-learner of lasso
#'
#' @description  X-learner, implemented via glmnet (lasso) with "coxph" distribution
#'
#' @param X The baseline covariates
#' @param W The treatment variable (0 or 1)
#' @param Y The follow-up time
#' @param D The event indicator
#' @param times The prediction time of interest
#' @param alpha Imbalance tuning parameter for a split (see grf documentation)
#' @param W.hat The propensity score
#' @param cen.fit The choice of model fitting for censoring
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' times <- 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[ ,1] + (-0.5 - 1 * X[ ,2]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.xl.lasso.fit <- surv_xl_lasso(X, Y, W, D, times, W.hat = 0.5)
#' cate <- predict(surv.xl.lasso.fit)
#' cate.test <- predict(surv.xl.lasso.fit, X.test)
#' }
#' @return A surv_xl_lasso object
#' @export
surv_xl_lasso <- function(X, Y, W, D, times, alpha = 0.05, W.hat = NULL, cen.fit = "Kaplan-Meier") {

  # fit model on W == 1 (cross-fitted using "preval" in glmnet)
  foldid1 <- sample(rep(seq(10), length = length(Y[W == 1])))
  x1 <- as.matrix(data.frame(X[W == 1, ]))
  lasso.fit1 <- glmnet::cv.glmnet(x1,
                                  survival::Surv(Y[W == 1], D[W == 1]),
                                  family = "cox",
                                  alpha = 1,
                                  keep = TRUE,
                                  foldid = foldid1)
  lambda.1.min <- lasso.fit1$lambda[which.min(lasso.fit1$cvm[!is.na(colSums(lasso.fit1$fit.preval))])]
  S0 <- base_surv(lasso.fit1, Y[W == 1], D[W == 1], x1, lambda = lambda.1.min)
  surf1 <- rep(NA, length(W))
  surf1[W == 1] <- pred_surv_preval(lasso.fit1, S0, times = times, lambda = lambda.1.min)
  surf1[W == 0] <- pred_surv(fit = lasso.fit1,
                             S0 = S0,
                             X = X[W == 0, ],
                             times = times,
                             lambda = lasso.fit1$lambda.min)

  # fit model on W == 0 (cross-fitted using "preval" in glmnet)
  foldid0 <- sample(rep(seq(10), length = length(Y[W == 0])))
  x0 <- as.matrix(data.frame(X[W == 0, ]))
  lasso.fit0 <- glmnet::cv.glmnet(x0,
                                  survival::Surv(Y[W == 0], D[W == 0]),
                                  family = "cox",
                                  alpha = 1,
                                  keep = TRUE,
                                  foldid = foldid0)
  lambda.0.min <- lasso.fit0$lambda[which.min(lasso.fit0$cvm[!is.na(colSums(lasso.fit0$fit.preval))])]
  S0 <- base_surv(lasso.fit0, Y[W == 0], D[W == 0], x0, lambda = lambda.0.min)
  surf0 <- rep(NA, length(W))
  surf0[W == 0] <- pred_surv_preval(lasso.fit0, S0, times = times, lambda = lambda.0.min)
  surf0[W == 1] <- pred_surv(fit = lasso.fit0,
                             S0 = S0,
                             X = X[W == 1, ],
                             times = times,
                             lambda = lasso.fit0$lambda.min)

  Tlasso1 <- 1 - surf1
  Tlasso0 <- 1 - surf0

  # IPCW weights (cross-fitted using "preval" in glmnet)
  if (cen.fit == "Kaplan-Meier") {
    shuffle <- sample(length(Y))
    kmdat <- data.frame(Y = Y[shuffle], D = D[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
    C.hat <- rep(NA, nrow(kmdat))
    for (z in 1:10) {
      testIndexes <- which(folds == z, arr.ind = TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      c.fit <- survival::survfit(survival::Surv(trainData$Y, 1 - trainData$D) ~ 1)
      cent <- testData$Y
      cent[testData$D == 0] <- times
      C.hat[testIndexes] <- summary(c.fit, times = cent)$surv
    }
    shudat <- data.frame(shuffle, C.hat)
    C.hat <- shudat[order(shuffle), ]$C.hat
  } else if (cen.fit == "survival.forest") {
    c.fit <- grf::survival_forest(cbind(W, X),
                                  Y,
                                  1 - D,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
    C.hat <- predict(c.fit)$predictions
    cent <- Y; cent[D == 0] <- times
    cen.times.index <- findInterval(cent, c.fit$failure.times)
    C.hat <- C.hat[cbind(1:length(Y), cen.times.index)]
  }

  # Propensity score (cross-fitted using "preval" in glmnet)
  if (is.null(W.hat)) {
    stop("propensity score needs to be supplied")
  } else {
    W.hat <- rep(W.hat, length(W))
  }

  sample.weights <- (1 / C.hat) * (1 / W.hat)

  # X-learner
  tempdat <- data.frame(Y = Y, D = D, W = W, sample.weights, X, Tlasso0, Tlasso1)
  binary.data <- tempdat[tempdat$D == 1 | tempdat$Y > times, ]
  binary.data$D[binary.data$D == 1 & binary.data$Y > times] <- 0
  binary.data <- binary.data[complete.cases(binary.data), ]
  b.data <- list(Y = binary.data$Y, D = binary.data$D, W = binary.data$W,
                 X = as.matrix(binary.data[ ,5:(ncol(binary.data)-2)]),
                 sample.weights = binary.data$sample.weights, mu0 = binary.data$Tlasso0, mu1 = binary.data$Tlasso1)

  foldid <- sample(rep(seq(10), length = length(b.data$Y[b.data$W == 1])))
  tau.fit1 <- glmnet::cv.glmnet(b.data$X[b.data$W == 1, ],
                                b.data$D[b.data$W == 1] - b.data$mu0[b.data$W == 1],
                                weights = b.data$sample.weights[b.data$W == 1],
                                foldid = foldid,
                                alpha = 1)
  XLtau1 <- as.vector(-predict(tau.fit1, X, s = "lambda.min"))

  foldid <- sample(rep(seq(10), length = length(b.data$Y[b.data$W == 0])))
  tau.fit0 <- glmnet::cv.glmnet(b.data$X[b.data$W == 0, ],
                                b.data$mu1[b.data$W == 0] - b.data$D[b.data$W == 0],
                                weights = b.data$sample.weights[b.data$W == 0],
                                foldid = foldid,
                                alpha = 1)
  XLtau0 <- as.vector(-predict(tau.fit0, X, s = "lambda.min"))

  # weighted CATE
  tau.hat <- as.vector(XLtau1 * (1 - W.hat) + XLtau0 * W.hat)

  ret <- list(tau.fit1 = tau.fit1,
              tau.fit0 = tau.fit0,
              W.hat = W.hat,
              tau.hat = tau.hat)
  class(ret) <- "surv_xl_lasso"
  ret
}

#' predict for surv_xl_lasso
#'
#' get estimated tau(X) using the trained surv_xl_lasso model
#'
#' @param object An surv_xl_lasso object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param W.hat The propensity score
#' @param ... Additional arguments (currently not used)
#'
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' times <- 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[ ,1] + (-0.5 - 1 * X[ ,2]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.xl.lasso.fit <- surv_xl_lasso(X, Y, W, D, times, W.hat = 0.5)
#' cate <- predict(surv.xl.lasso.fit)
#' cate.test <- predict(surv.xl.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_xl_lasso <- function(object,
                                  newdata = NULL,
                                  W.hat = NULL,
                                  ...) {
  if (is.null(newdata)) {
    return(object$tau.hat)
  } else {
    XLtau1 <- as.vector(-predict(object$tau.fit1, newdata, s = "lambda.min"))
    XLtau0 <- as.vector(-predict(object$tau.fit0, newdata, s = "lambda.min"))
    if (is.null(W.hat)) {
      W.hat <- rep(object$W.hat[1], nrow(newdata))
    } else {
      if (length(W.hat) == 1) {
      W.hat <- rep(W.hat, nrow(newdata))
      }
    }
    return(as.vector(XLtau1 * (1 - W.hat) + XLtau0 * W.hat))
  }
}
