#' @title X-learner of grf
#'
#' @description  X-learner, implemented via survival_forest in the grf package
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
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
#' surv.xl.grf.fit <- surv_xl_grf(X, Y, W, D, times, W.hat = 0.5)
#' cate <- predict(surv.xl.grf.fit)
#' cate.test <- predict(surv.xl.grf.fit, X.test)
#' }
#' @return A surv_xl_grf object
#' @export
surv_xl_grf <- function(X, Y, W, D, times, alpha = 0.05, W.hat = NULL, cen.fit = "Kaplan-Meier") {
  # fit model on W == 1
  grffit1 <- grf::survival_forest(X[W == 1, ],
                                  Y[W == 1],
                                  D[W == 1],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
  surf1 <- rep(NA, length(W))
  times.index <- findInterval(times, grffit1$failure.times)
  surf1[W == 1] <- predict(grffit1)$predictions[ ,times.index]
  surf1[W == 0] <- predict(grffit1, X[W == 0, ])$predictions[ ,times.index]

  # fit model on W == 0
  grffit0 <- grf::survival_forest(X[W == 0, ],
                                  Y[W == 0],
                                  D[W == 0],
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
  surf0 <- rep(NA, length(W))
  times.index <- findInterval(times, grffit0$failure.times)
  surf0[W == 0] <- predict(grffit0)$predictions[ ,times.index]
  surf0[W == 1] <- predict(grffit0, X[W == 1, ])$predictions[ ,times.index]

  Tgrf1 <- 1 - surf1
  Tgrf0 <- 1 - surf0

  # IPCW weights
  Q <- as.numeric(D == 1 | Y > times)    # indicator for uncensored at t0
  U <- pmin(Y, times)                         # truncated follow-up time by t0
  if (cen.fit == "Kaplan-Meier") {
    shuffle <- sample(length(U))
    kmdat <- data.frame(U = U[shuffle], Q = Q[shuffle])
    folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
    C.hat <- rep(NA, nrow(kmdat))
    for (z in 1:10) {
      testIndexes <- which(folds == z, arr.ind = TRUE)
      testData <- kmdat[testIndexes, ]
      trainData <- kmdat[-testIndexes, ]
      c.fit <- survival::survfit(survival::Surv(trainData$U, 1 - trainData$Q) ~ 1)
      C.hat[testIndexes] <- summary(c.fit, times = testData$U)$surv
    }
    shudat <- data.frame(shuffle, C.hat)
    C.hat <- shudat[order(shuffle), ]$C.hat
  } else if (cen.fit == "survival.forest") {
    c.fit <- grf::survival_forest(cbind(W, X),
                                  U,
                                  1 - Q,
                                  alpha = alpha,
                                  prediction.type = "Nelson-Aalen")
    C.hat <- predict(c.fit)$predictions
    cen.times.index <- findInterval(U, c.fit$failure.times)
    C.hat <- C.hat[cbind(1:length(U), cen.times.index)]
  }
  ipcw <- 1 / C.hat

  # Propensity score
  if (is.null(W.hat)) {
    stop("propensity score needs to be supplied")
  } else {
    W.hat <- rep(W.hat, length(W))
  }

  sample.weights <- ipcw / W.hat  # censoring weight * treatment weight

  # X-learner
  tempdat <- data.frame(Y = Y, D = D, W = W, sample.weights, X, Tgrf0, Tgrf1)
  binary.data <- tempdat[tempdat$D == 1 | tempdat$Y > times, ]
  binary.data$D[binary.data$D == 1 & binary.data$Y > times] <- 0
  binary.data <- binary.data[complete.cases(binary.data), ]
  b.data <- list(Y = binary.data$Y, D = binary.data$D, W = binary.data$W,
                 X = as.matrix(binary.data[ ,5:(ncol(binary.data)-2)]),
                 sample.weights = binary.data$sample.weights, mu0 = binary.data$Tgrf0, mu1 = binary.data$Tgrf1)

  tau.fit1 <- grf::regression_forest(b.data$X[b.data$W == 1, ],
                                     b.data$D[b.data$W == 1] - b.data$mu0[b.data$W == 1],
                                     sample.weights = b.data$sample.weights[b.data$W == 1])
  XLtau1 <- -predict(tau.fit1, data.frame(X))

  tau.fit0 <- grf::regression_forest(b.data$X[b.data$W == 0, ],
                                     b.data$mu1[b.data$W == 0] - b.data$D[b.data$W == 0],
                                     sample.weights = b.data$sample.weights[b.data$W == 0])
  XLtau0 <- -predict(tau.fit0, data.frame(X))

  # weighted CATE
  tau.hat <- as.vector(XLtau1 * (1 - W.hat) + XLtau0 * W.hat)

  ret <- list(tau.fit1 = tau.fit1,
              tau.fit0 = tau.fit0,
              W.hat = W.hat,
              tau.hat = tau.hat)
  class(ret) <- "surv_xl_grf"
  ret
}

#' predict for surv_xl_grf
#'
#' get estimated tau(X) using the trained surv_xl_grf model
#'
#' @param object An surv_xl_grf object
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
#' surv.xl.grf.fit <- surv_xl_grf(X, Y, W, D, times, W.hat = 0.5)
#' cate <- predict(surv.xl.grf.fit)
#' cate.test <- predict(surv.xl.grf.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_xl_grf <- function(object,
                                newdata = NULL,
                                W.hat = NULL,
                                ...) {
  if (is.null(newdata)) {
    return(object$tau.hat)
  } else {
    XLtau1 <- -predict(object$tau.fit1, data.frame(newdata))
    XLtau0 <- -predict(object$tau.fit0, data.frame(newdata))
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
