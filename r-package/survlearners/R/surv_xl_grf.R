#' @title X-learner of grf
#'
#' @description  X-learner, implemented via survival_forest in the grf package
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param W.hat The propensity score
#' @param cen.fit The choice of model fitting for censoring
#' @param k.folds The number of folds for estimating nuisance parameters via cross-fitting
#' @param args.grf.nuisance Input arguments for a grf model that estimates nuisance parameters
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' t0 <- 0.2
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
#' surv.xl.grf.fit <- surv_xl_grf(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.xl.grf.fit)
#' cate.test <- predict(surv.xl.grf.fit, X.test)
#' }
#' @return A surv_xl_grf object
#' @export
surv_xl_grf <- function(X, Y, W, D, t0, W.hat = NULL, cen.fit = "Kaplan-Meier",
                        k.folds = 10, args.grf.nuisance = list()) {

  args.grf.nuisance <- list(failure.times = NULL,
                            num.trees = max(50, 2000 / 4),
                            sample.weights = NULL,
                            clusters = NULL,
                            equalize.cluster.weights = FALSE,
                            sample.fraction = 0.5,
                            mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                            min.node.size = 15,
                            honesty = TRUE,
                            honesty.fraction = 0.5,
                            honesty.prune.leaves = TRUE,
                            alpha = 0.05,
                            prediction.type = "Nelson-Aalen",
                            compute.oob.predictions = TRUE,
                            num.threads = NULL,
                            seed = runif(1, 0, .Machine$integer.max))

  # fit model on W == 1
  grffit1 <- do.call(grf::survival_forest, c(list(X = X[W == 1, ], Y = Y[W == 1], D = D[W == 1]), args.grf.nuisance))
  surf1 <- rep(NA, length(W))
  t0.index <- findInterval(t0, grffit1$failure.times)
  if (t0.index == 0) {
    surf1[W == 1] <- rep(1, length(D[W == 1]))
    surf1[W == 0] <- rep(1, length(D[W == 0]))
  } else {
    surf1[W == 1] <- predict(grffit1)$predictions[ ,t0.index]
    surf1[W == 0] <- predict(grffit1, X[W == 0, ])$predictions[ ,t0.index]
  }

  # fit model on W == 0
  grffit0 <- do.call(grf::survival_forest, c(list(X = X[W == 0, ], Y = Y[W == 0], D = D[W == 0]), args.grf.nuisance))
  surf0 <- rep(NA, length(W))
  t0.index <- findInterval(t0, grffit0$failure.times)
  if (t0.index == 0) {
    surf0[W == 0] <- rep(1, length(D[W == 0]))
    surf0[W == 1] <- rep(1, length(D[W == 1]))
  } else {
    surf0[W == 0] <- predict(grffit0)$predictions[ ,t0.index]
    surf0[W == 1] <- predict(grffit0, X[W == 1, ])$predictions[ ,t0.index]
  }

  Tgrf1 <- 1 - surf1
  Tgrf0 <- 1 - surf0

  # IPCW weights
  Q <- as.numeric(D == 1 | Y > t0)    # indicator for uncensored at t0
  U <- pmin(Y, t0)                         # truncated follow-up time by t0
  if (cen.fit == "Kaplan-Meier") {
    fold.id <- sample(rep(seq(k.folds), length = nrow(X)))
    C.hat <- rep(NA, length(fold.id))
    for (z in 1:k.folds) {
      c.fit <- survival::survfit(survival::Surv(Y[!fold.id == z], 1 - Q[!fold.id == z]) ~ 1)
      C.hat[fold.id == z] <- summary(c.fit, times = U[fold.id == z])$surv
    }
  } else if (cen.fit == "survival.forest") {
    c.fit <- do.call(grf::survival_forest, c(list(X = cbind(W, X), Y = Y, D = 1 - Q), args.grf.nuisance))
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
  binary.data <- tempdat[tempdat$D == 1 | tempdat$Y > t0, ]
  binary.data$D[binary.data$D == 1 & binary.data$Y > t0] <- 0
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
#' t0 <- 0.2
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
#' surv.xl.grf.fit <- surv_xl_grf(X, Y, W, D, t0, W.hat = 0.5)
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
