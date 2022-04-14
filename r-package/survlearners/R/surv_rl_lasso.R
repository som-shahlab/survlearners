#' @title R-learner of lasso
#'
#' @description  R-learner, implemented via glmnet (lasso)
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param times The prediction time of interest
#' @param k.folds Number of folds for cross validation
#' @param foldid User-supplied foldid. Must have length equal to length(W). If provided, it overrides the k.folds option.
#' @param W.hat Propensity score
#' @param Y.hat Conditional mean outcome E(Y|X)
#' @param C.hat Censoring weights
#' @param lambda.choice How to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty.factor User-supplied penalty factor, must be of length the same as number of features in X
#' @param args.lasso.nuisance Input arguments for a lasso model that estimates nuisance parameters
#' @param args.grf.nuisance Input arguments for a grf model that estimates nuisance parameters
#' @param args.lasso.tau Input arguments for a lasso model that estimates CATE
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
#' surv.rl.lasso.fit <- surv_rl_lasso(X, Y, W, D, times, W.hat = 0.5)
#' cate <- predict(surv.rl.lasso.fit)
#' cate.test <- predict(surv.rl.lasso.fit, X.test)
#' }
#' @return A surv_rl_lasso object
#' @export
surv_rl_lasso <- function(X, Y, W, D,
                          times = NULL,
                          k.folds = 10,
                          foldid = NULL,
                          W.hat = NULL,
                          Y.hat = NULL,
                          C.hat = NULL,
                          lambda.choice = "lambda.min",
                          penalty.factor = NULL,
                          args.lasso.nuisance = list(),
                          args.grf.nuisance = list(),
                          args.lasso.tau = list(),
                          cen.fit = "Kaplan-Meier") {

    input <- sanitize_input(X, Y, W, D)
    X <- input$X
    W <- input$W
    Y <- input$Y
    D <- input$D

    x.scl <- scale(X, center = TRUE, scale = TRUE)
    x.scl <- x.scl[ ,!is.na(colSums(x.scl)), drop = FALSE]

    nobs <- nrow(x.scl)
    pobs <- ncol(x.scl)

    if (is.null(foldid) || length(foldid) != length(W)) {

      if (!is.null(foldid) && length(foldid) != length(W)) {
        warning("supplied foldid does not have the same length ")
      }

      if (is.null(k.folds)) {
          k.folds <- floor(max(3, min(10, length(W) / 4)))
      }

      # fold ID for cross-validation; balance treatment assignments
      foldid <- sample(rep(seq(k.folds), length = length(W)))

    }

    # penalty factor for nuisance and tau estimators
    if (is.null(penalty.factor) || (length(penalty.factor) != pobs)) {
      if (!is.null(penalty.factor) && length(penalty.factor) != pobs) {
        warning("penalty.factor supplied is not of the same length as the number of columns in X after removing NA columns. Using all ones instead.")
      }
      penalty.factor.nuisance.w <- rep(1, pobs)
      penalty.factor.nuisance.m <- rep(1, (pobs+1))
      penalty.factor.tau <- c(0, rep(1, pobs))
    } else {
      penalty.factor.nuisance <- penalty.factor
      penalty.factor.tau <- c(0, penalty.factor)
    }

    if (is.null(W.hat)) {
      stop("propensity score needs to be supplied")
    } else if (length(W.hat) == 1) {
      W.hat <- rep(W.hat, nrow(X))
    } else if (length(W.hat) != nrow(X)) {
      stop("W.hat has incorrect length.")
    }


    args.lasso.nuisance <- list(family = "cox",
                                nfolds = 10,
                                alpha = 1,
                                lambda = NULL,
                                penalty.factor = penalty.factor.nuisance.m)

    if (is.null(Y.hat)) {
    foldid <- sample(rep(seq(k.folds), length = length(W)))
    survt1 <- survt0 <- rep(NA, length(W))
    for (k in 1:k.folds) {
      XW <- as.matrix(data.frame(W[!foldid == k], X[!foldid == k, ]))
      y <- survival::Surv(Y[!foldid == k], D[!foldid == k])
      y.fit <- do.call(glmnet::cv.glmnet, c(list(x = XW, y = y), args.lasso.nuisance))
      S0 <- base_surv(y.fit, Y[!foldid == k], D[!foldid == k], XW, lambda = y.fit$lambda.min)
      survt1[foldid == k] <- pred_surv(y.fit, S0, cbind(rep(1, length(W[foldid == k])), X[foldid == k, ]), times = times, lambda = y.fit$lambda.min)
      survt0[foldid == k] <- pred_surv(y.fit, S0, cbind(rep(0, length(W[foldid == k])), X[foldid == k, ]), times = times, lambda = y.fit$lambda.min)
    }
    Y.hat  <- W.hat * survt1 + (1 - W.hat) * survt0
    } else {
      y.fit <- NULL
    }

    args.grf.nuisance <- list(failure.times = NULL,
                              num.trees = max(50, 2000 / 4),
                              min.node.size = 15,
                              honesty = TRUE,
                              honesty.fraction = 0.5,
                              honesty.prune.leaves = TRUE,
                              alpha = 0.05,
                              prediction.type = "Nelson-Aalen",
                              compute.oob.predictions = TRUE)

    if (is.null(C.hat)) {
      Q <- as.numeric(D == 1 | Y > times)         # indicator for uncensored at t0
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
    } else {
      c.fit <- NULL
    }

    # use binary data
    tempdat <- data.frame(Y, D, W, Y.hat, W.hat, C.hat, foldid, x.scl)
    binary.data <- tempdat[tempdat$D == 1 | tempdat$Y > times, ]          # remove subjects who got censored before the time of interest t50
    binary.data$D[binary.data$D == 1 & binary.data$Y > times] <- 0     # recode the event status for subjects who had events after t50
    binary.data <- binary.data[complete.cases(binary.data), ]

    sample.weights <- 1 / binary.data$C.hat
    y.tilde <- (1 - binary.data$D) - binary.data$Y.hat
    x.scl <- binary.data[ , 8:dim(binary.data)[2]]
    foldid2 <- sample(rep(seq(k.folds), length = length(binary.data$W)))

    x.scl.tilde <- cbind(as.numeric(binary.data$W - binary.data$W.hat) * cbind(1, x.scl))
    x.scl.pred <- cbind(1, x.scl)

    args.lasso.tau <- list(weights = sample.weights,
                           foldid = foldid2,
                           alpha = 1,
                           lambda = NULL,
                           penalty.factor = penalty.factor.tau,
                           standardize = FALSE)
    tau.fit <- do.call(glmnet::cv.glmnet, c(list(x = as.matrix(x.scl.tilde), y = y.tilde), args.lasso.tau))
    tau.beta <- as.vector(t(coef(tau.fit, s = lambda.choice)[-1]))
    tau.hat <- as.matrix(x.scl.pred) %*% tau.beta

    ret <- list(tau.fit = tau.fit,
               tau.beta = tau.beta,
               y.fit = y.fit,
               c.fit = c.fit,
               W.hat = W.hat,
               Y.hat = Y.hat,
               C.hat = C.hat,
               tau.hat = tau.hat)
    class(ret) <- "surv_rl_lasso"
    ret
}


#' predict for surv_rl_lasso
#'
#' get estimated tau(X) using the trained surv_rl_lasso model
#'
#' @param object An surv_rl_lasso object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
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
#' surv.rl.lasso.fit <- surv_rl_lasso(X, Y, W, D, times, W.hat = 0.5)
#' cate <- predict(surv.rl.lasso.fit)
#' cate.test <- predict(surv.rl.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_rl_lasso <- function(object,
                                  newdata = NULL,
                                  ...) {
  if (!is.null(newdata)) {
    newdata <- sanitize_x(newdata)
    newdata.scl <- scale(newdata, center = TRUE, scale = TRUE)
    newdata.scl <- newdata.scl[ ,!is.na(colSums(newdata.scl)), drop = FALSE]
    newdata.scl.pred <- cbind(1, newdata.scl)
    tau.hat <- newdata.scl.pred %*% object$tau.beta
  } else {
    tau.hat <- object$tau.hat
  }
  return(tau.hat)
}
