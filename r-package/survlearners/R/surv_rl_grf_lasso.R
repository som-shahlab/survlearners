#' @title R-learner of grf and lasso
#'
#' @description  R-learner, implemented via the grf package for nuisance parameter estimation and lasso for target parameter
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param times The prediction time of interest
#' @param k.folds Number of folds for cross validation
#' @param W.hat Propensity score
#' @param Y.hat Conditional mean outcome E(Y|X)
#' @param C.hat Censoring weights
#' @param penalty.factor User-supplied penalty factor, must be of length the same as number of features in X
#' @param lambda.choice How to cross-validate; choose from "lambda.min" or "lambda.1se"
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
#' surv.rl.grf.lasso.fit <- surv_rl_grf_lasso(X, Y, W, D, times, W.hat = 0.5)
#' cate <- predict(surv.rl.grf.lasso.fit)
#' cate.test <- predict(surv.rl.grf.lasso.fit, X.test)
#' }
#' @return A surv_rl_grf_lasso object
#' @export
surv_rl_grf_lasso <- function(X, Y, W, D,
                              times = NULL,
                              k.folds = 10,
                              W.hat = NULL,
                              Y.hat = NULL,
                              C.hat = NULL,
                              penalty.factor = NULL,
                              lambda.choice = "lambda.min",
                              args.grf.nuisance = list(),
                              args.lasso.tau = list(),
                              cen.fit = "Kaplan-Meier") {

  input <- sanitize_input(X, Y, W, D)
  X <- input$X
  W <- as.numeric(input$W)
  Y <- input$Y
  D <- input$D
  nobs <- nrow(X)
  pobs <- ncol(X)

  x.scl <- scale(X, center = TRUE, scale = TRUE)
  x.scl <- x.scl[ ,!is.na(colSums(x.scl)), drop = FALSE]

  # penalty factor for tau estimator
  if (is.null(penalty.factor) || (length(penalty.factor) != pobs)) {
    if (!is.null(penalty.factor) && length(penalty.factor) != pobs) {
      warning("penalty.factor supplied is not of the same length as the number of columns in X after removing NA columns. Using all ones instead.")
    }
    penalty.factor.tau <- c(0, rep(1, pobs))
  } else {
    penalty.factor.tau <- c(0, penalty.factor)
  }

  if (is.null(W.hat)) {
    stop("propensity score needs to be supplied")
  } else if (length(W.hat) == 1) {
    W.hat <- rep(W.hat, nrow(X))
  } else if (length(W.hat) != nrow(X)) {
    stop("W.hat has incorrect length.")
  }

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

  if (is.null(Y.hat)) {
    y.fit <- do.call(grf::survival_forest, c(list(X = cbind(X, W), Y = Y, D = D), args.grf.nuisance))
    y.fit[["X.orig"]][ , ncol(X) + 1] <- rep(1, nrow(X))
    S1.hat <- predict(y.fit)$predictions
    y.fit[["X.orig"]][ , ncol(X) + 1] <- rep(0, nrow(X))
    S0.hat <- predict(y.fit)$predictions
    y.fit[["X.orig"]][ , ncol(X) + 1] <- W

    times.index <- findInterval(times, y.fit$failure.times)
    surf1 <- S1.hat[ , times.index]
    surf0 <- S0.hat[ , times.index]
    Y.hat  <- W.hat * surf1 + (1 - W.hat) * surf0
  } else {
    y.fit <- NULL
  }


  if (is.null(C.hat)) {
    if (cen.fit == "Kaplan-Meier") {
      traindat <- data.frame(Y = Y, D = D)
      shuffle <- sample(nrow(traindat))
      kmdat <- traindat[shuffle, ]
      folds <- cut(seq(1, nrow(kmdat)), breaks = 10, labels = FALSE)
      C.hat <- rep(NA, nrow(kmdat))
      for (z in 1:10) {
        testIndexes <- which(folds == z, arr.ind = TRUE)
        testData <- kmdat[testIndexes, ]
        trainData <- kmdat[-testIndexes, ]
        c.fit <- survival::survfit(survival::Surv(trainData$Y, 1 - trainData$D) ~ 1)
        cent <- testData$Y; cent[testData$D == 0] <- times
        C.hat[testIndexes] <- summary(c.fit, times = cent)$surv
      }
      shudat <- data.frame(shuffle, C.hat)
      C.hat <- shudat[order(shuffle), ]$C.hat
    } else if (cen.fit == "survival.forest") {
      c.fit <- do.call(grf::survival_forest, c(list(X = cbind(X, W), Y = Y, D = 1 - D), args.grf.nuisance))
      C.hat <- predict(c.fit, failure.times = c.fit$failure.times)$predictions
      cent <- Y; cent[D == 0] <- times
      cen.times.index <- findInterval(cent, c.fit$failure.times)
      C.hat <- C.hat[cbind(1:length(Y), cen.times.index)]
    }
  } else {
    c.fit <- NULL
  }

  # create binary data
  tempdat <- data.frame(Y, D, W, Y.hat, W.hat, C.hat, x.scl)
  binary.data <- tempdat[tempdat$D == 1 | tempdat$Y > times, ]          # remove subjects who got censored before the time of interest t50
  binary.data$D[binary.data$D == 1 & binary.data$Y > times] <- 0     # recode the event status for subjects who had events after t50
  binary.data <- binary.data[complete.cases(binary.data), ]

  sample.weights <- 1 / binary.data$C.hat     # the treatment weight is already accounted
  y.tilde <- (1 - binary.data$D) - binary.data$Y.hat
  x.scl <- binary.data[ , 7:dim(binary.data)[2]]
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
              tau.hat = tau.hat)
  class(ret) <- "surv_rl_grf_lasso"
  ret
}

#' predict for surv_rl_grf_lasso
#'
#' get estimated tau(X) using the trained surv_rl_grf_lasso model
#'
#' @param object A surv_rl_grf_lasso object
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
#' surv.rl.grf.lasso.fit <- surv_rl_grf_lasso(X, Y, W, D, times, W.hat = 0.5)
#' cate <- predict(surv.rl.grf.lasso.fit)
#' cate.test <- predict(surv.rl.grf.lasso.fit, X.test)
#' }
#'
#' @return A vector of predicted conditional average treatment effects
#' @export
predict.surv_rl_grf_lasso <- function(object,
                                      newdata = NULL,
                                      ...) {
  if (!is.null(newdata)) {
    newdata <- sanitize_x(newdata)
    newdata.scl <- scale(newdata, center = TRUE, scale = TRUE)
    newdata.scl <- newdata.scl[ ,!is.na(colSums(newdata.scl)), drop = FALSE]
    newdata.scl.pred <- cbind(1, newdata.scl)
    tau.hat <- newdata.scl.pred %*% object$tau.beta
  }
  else {
    tau.hat <- object$tau.hat
  }
  return(tau.hat)
}
