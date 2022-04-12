#' @title R-learner of grf
#'
#' @description  R-learner, implemented via the grf package
#'
#' @param X The baseline covariates
#' @param W The treatment variable (0 or 1)
#' @param Y The follow-up time
#' @param D The event indicator
#' @param times The prediction time of interest
#' @param k.folds Number of folds for cross validation
#' @param W.hat Propensity score
#' @param Y.hat Conditional mean outcome E(Y|X)
#' @param C.hat Censoring weights
#' @param args.grf.nuisance Input arguments for a grf model that estimates nuisance parameters
#' @param args.grf.tau Input arguments for a grf model that estimates CATE
#' @param cen.fit The choice of model fitting for censoring
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
#' surv.rl.grf.fit = surv_rl_grf(X, Y, W, D, times, W.hat = 0.5)
#' cate = predict(surv.rl.grf.fit)
#' cate.test = predict(surv.rl.grf.fit, X.test)
#' }
#' @return a surv_rl_grf_fit object
#' @export
surv_rl_grf = function(X, Y, W, D,
                       times = NULL,
                       k.folds = NULL,
                       W.hat = NULL,
                       Y.hat = NULL,
                       C.hat = NULL,
                       args.grf.nuisance = list(),
                       args.grf.tau = list(),
                       cen.fit = "Kaplan-Meier"){

  input = sanitize_input(X, Y, W, D)
  X = input$X
  W = as.numeric(input$W)
  Y = input$Y
  D = input$D
  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(k.folds)) {
    k.folds = floor(max(3, min(10,length(Y)/4)))
  }

  if (is.null(W.hat)){
    stop("propensity score needs to be supplied")
  }else if (length(W.hat) == 1) {
    W.hat <- rep(W.hat, nrow(X))
  }else if (length(W.hat) != nrow(X)){
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

  if (is.null(Y.hat)){
    y.fit <- do.call(grf::survival_forest, c(list(X = cbind(X, W), Y = Y, D = D), args.grf.nuisance))
    y.fit[["X.orig"]][, ncol(X) + 1] <- rep(1, nrow(X))
    S1.hat <- predict(y.fit)$predictions
    y.fit[["X.orig"]][, ncol(X) + 1] <- rep(0, nrow(X))
    S0.hat <- predict(y.fit)$predictions
    y.fit[["X.orig"]][, ncol(X) + 1] <- W

    times.index <- findInterval(times, y.fit$failure.times)
    surf1 <- S1.hat[, times.index]
    surf0 <- S0.hat[, times.index]
    Y.hat  <- W.hat * surf1 + (1 - W.hat) * surf0
  }else {
    y.fit = NULL
  }

  if (is.null(C.hat)){
    if(cen.fit == "Kaplan-Meier"){
      traindat <- data.frame(Y = Y, D = D)
      shuffle <- sample(nrow(traindat))
      kmdat <- traindat[shuffle,]
      folds <- cut(seq(1, nrow(kmdat)), breaks=10, labels=FALSE)
      C.hat <- rep(NA, nrow(kmdat))
      for(z in 1:10){
        testIndexes <- which(folds==z, arr.ind=TRUE)
        testData <- kmdat[testIndexes, ]
        trainData <- kmdat[-testIndexes, ]
        c.fit <- survival::survfit(survival::Surv(trainData$Y, 1 - trainData$D) ~ 1)
        cent <- testData$Y; cent[testData$D==0] <- times
        C.hat[testIndexes] <- summary(c.fit, times = cent)$surv
      }
      shudat <- data.frame(shuffle, C.hat)
      C.hat <- shudat[order(shuffle), ]$C.hat
    }else if (cen.fit == "survival.forest"){
      c.fit <- do.call(grf::survival_forest, c(list(X = cbind(X, W), Y = Y, D = 1 - D), args.grf.nuisance))
      C.hat <- predict(c.fit, failure.times = c.fit$failure.times)$predictions
      cent <- Y; cent[D==0] <- times
      cen.times.index <- findInterval(cent, c.fit$failure.times)
      C.hat <- C.hat[cbind(1:length(Y), cen.times.index)]
    }
  }else{
    c.fit <- NULL
  }

  # create binary data
  tempdata <- data.frame(Y, D, W, Y.hat, W.hat, C.hat, X)
  binary.data <- tempdata[tempdata$D==1|tempdata$Y > times,]       # remove subjects who got censored before the time of interest t50
  binary.data$D[binary.data$D==1 & binary.data$Y > times] <- 0     # recode the event status for subjects who had events after t50
  binary.data <- binary.data[complete.cases(binary.data),]

  y.tilde <- (1 - binary.data$D) - binary.data$Y.hat
  w.tilde <-  binary.data$W - binary.data$W.hat
  pseudo.outcome <- y.tilde/w.tilde
  sample.weights <- w.tilde^2/binary.data$C.hat

  args.grf.tau <- list(sample.weights = sample.weights,
                            num.trees = 2000,
                            clusters = NULL,
                            sample.fraction = 0.5,
                            mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                            min.node.size = 5,
                            honesty = TRUE,
                            honesty.fraction = 0.5,
                            honesty.prune.leaves = TRUE,
                            imbalance.penalty = 0,
                            ci.group.size = 2,
                            compute.oob.predictions = TRUE,
                            num.threads = NULL,
                            seed = runif(1, 0, .Machine$integer.max))

  tau.fit <- do.call(grf::regression_forest, c(list(X = binary.data[,7:dim(binary.data)[2]], Y = pseudo.outcome), args.grf.tau))

  ret <- list(tau.fit = tau.fit,
              pseudo.outcome = pseudo.outcome,
              sample.weights = sample.weights,
              y.fit = y.fit,
              c.fit = c.fit,
              W.hat = W.hat,
              Y.hat = Y.hat,
              C.hat = C.hat)
  class(ret) <- "surv_rl_grf"
  ret
}

#' predict for surv_rl_grf
#'
#' get estimated tau(X) using the trained surv_rl_grf model
#'
#' @param object a surv_rl_grf object
#' @param newdata covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param tau.only if set to TRUE, onlly return prediction on tau. Otherwise, return a list including prediction on tau, propensity score, and baseline main effect.
#' @param ... additional arguments (currently not used)
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
#' surv.rl.grf.fit = surv_rl_grf(X, Y, W, D, times, W.hat = 0.5)
#' cate = predict(surv.rl.grf.fit)
#' cate.test = predict(surv.rl.grf.fit, X.test)
#' }
#'
#' @return A vector of predicted conditional average treatment effects
#' @export
predict.surv_rl_grf <- function(object,
                                newdata = NULL,
                                tau.only = TRUE,
                                ...) {
  if (!is.null(newdata)){
    newdata = sanitize_x(newdata)
  }
  if (tau.only) {
    return(predict(object$tau.fit, newdata)$predictions)
  } else {
    tau <- predict(object$tau.fit, newdata)$predictions
    e = predict(object$w.fit, newdata)$predictions
    m = predict(object$y.fit, newdata)$predictions
    mu1 = m + (1-e) * tau
    mu0 = m - e * tau
    return(list(tau=tau, e=e, m=m, mu1 = mu1, mu0 = mu0))
  }
}
