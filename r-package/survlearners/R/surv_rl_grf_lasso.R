#' @title R-learner of grf and lasso
#'
#' @description  R-learner, implemented via the grf package for nuisance parameter estimation and lasso for target parameter
#'
#' @param X The baseline covariates
#' @param W The treatment variable (0 or 1)
#' @param Y The follow-up time
#' @param D The event indicator
#' @param k.folds Number of folds for cross validation
#' @param p.hat Propensity score
#' @param m.hat Conditional mean outcome E(Y|X)
#' @param c.hat Censoring weights
#' @param times The prediction time of interest
#' @param failure.times A vector of event times to fit the survival curve at.
#' @param num.trees Number of trees grown in the forest
#' @param alpha Imbalance tuning parameter for a split in a tree
#' @param sample.weights See grf documentation
#' @param clusters See grf documentation
#' @param equalize.cluster.weights See grf documentation
#' @param sample.fraction See grf documentation
#' @param mtry See grf documentation
#' @param min.node.size See grf documentation
#' @param honesty See grf documentation
#' @param honesty.fraction See grf documentation
#' @param honesty.prune.leaves See grf documentation
#' @param alpha See grf documentation
#' @param imbalance.penalty See grf documentation
#' @param stabilize.splits See grf documentation
#' @param ci.group.size See grf documentation
#' @param tune.parameters See grf documentation
#' @param compute.oob.predictions See grf documentation
#' @param num.threads See grf documentation
#' @param seed See grf documentation
#' @param lambda.tau User-supplied lambda sequence for cross validation in the cate model
#' @param lambda.choice How to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty.factor User-supplied penalty factor, must be of length the same as number of features in X
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
#' surv.rl.grf.lasso.fit = surv_rl_grf_lasso(X, W, Y, D, times, p.hat = 0.5)
#' cate = predict(surv.rl.grf.lasso.fit)
#' cate.test = predict(surv.rl.grf.lasso.fit, X.test)
#' }
#' @return a surv_rl_grf_lasso object
#' @export
surv_rl_grf_lasso = function(X, W, Y, D,
                             times = NULL,
                             k.folds = 10,
                             p.hat = NULL,
                             m.hat = NULL,
                             c.hat = NULL,
                             failure.times = NULL,
                             num.trees = 2000,
                             sample.weights = NULL,
                             clusters = NULL,
                             equalize.cluster.weights = FALSE,
                             sample.fraction = 0.5,
                             mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                             min.node.size = 5,
                             honesty = TRUE,
                             honesty.fraction = 0.5,
                             honesty.prune.leaves = TRUE,
                             alpha = 0.05,
                             imbalance.penalty = 0,
                             stabilize.splits = TRUE,
                             ci.group.size = 2,
                             tune.parameters = "none",
                             compute.oob.predictions = TRUE,
                             num.threads = NULL,
                             seed = runif(1, 0, .Machine$integer.max),
                             lambda.tau = NULL,
                             lambda.choice = "lambda.min",
                             penalty.factor = NULL,
                             cen.fit = "KM",
                             verbose = FALSE){

  input = sanitize_input(X,W,Y,D)
  X = input$X
  W = as.numeric(input$W)
  Y = input$Y
  D = input$D
  nobs = nrow(X)
  pobs = ncol(X)

  x.scl = scale(X, center = TRUE, scale = TRUE)
  x.scl = x.scl[,!is.na(colSums(x.scl)), drop = FALSE]

  # penalty factor for tau estimator
  if (is.null(penalty.factor) || (length(penalty.factor) != pobs)) {
    if (!is.null(penalty.factor) && length(penalty.factor) != pobs) {
      warning("penalty.factor supplied is not of the same length as the number of columns in X after removing NA columns. Using all ones instead.")
    }
    penalty.factor.tau = c(0, rep(1, pobs))
  } else {
    penalty.factor.tau = c(0, penalty.factor)
  }

  if (is.null(p.hat)){
    stop("propensity score needs to be supplied")
  }else if (length(p.hat) == 1) {
    p.hat <- rep(p.hat, nrow(X))
  }else if (length(p.hat) != nrow(X)){
    stop("p.hat has incorrect length.")
  }

  args.nuisance <- list(failure.times = failure.times,
                        num.trees = max(50, num.trees / 4),
                        sample.weights = sample.weights,
                        clusters = clusters,
                        equalize.cluster.weights = equalize.cluster.weights,
                        sample.fraction = sample.fraction,
                        mtry = mtry,
                        min.node.size = 15,
                        honesty = TRUE,
                        honesty.fraction = 0.5,
                        honesty.prune.leaves = TRUE,
                        alpha = alpha,
                        prediction.type = "Nelson-Aalen",
                        compute.oob.predictions = FALSE,
                        num.threads = num.threads,
                        seed = seed)

  if (is.null(m.hat)){
    y.fit <- do.call(grf::survival_forest, c(list(X = cbind(X, W), Y = Y, D = D), args.nuisance))
    y.fit[["X.orig"]][, ncol(X) + 1] <- rep(1, nrow(X))
    S1.hat <- predict(y.fit)$predictions
    y.fit[["X.orig"]][, ncol(X) + 1] <- rep(0, nrow(X))
    S0.hat <- predict(y.fit)$predictions
    y.fit[["X.orig"]][, ncol(X) + 1] <- W

    times.index <- findInterval(times, y.fit$failure.times)
    surf1 <- S1.hat[, times.index]
    surf0 <- S0.hat[, times.index]
    m.hat  <- p.hat * surf1 + (1 - p.hat) * surf0
  }else {
    y.fit = NULL
  }

  args.nuisance$compute.oob.predictions <- TRUE
  if (is.null(c.hat)){
    if(cen.fit == "KM"){
      traindat <- data.frame(Y = Y, D = D)
      shuffle <- sample(nrow(traindat))
      kmdat <- traindat[shuffle,]
      folds <- cut(seq(1, nrow(kmdat)), breaks=10, labels=FALSE)
      c.hat <- rep(NA, nrow(kmdat))
      for(z in 1:10){
        testIndexes <- which(folds==z, arr.ind=TRUE)
        testData <- kmdat[testIndexes, ]
        trainData <- kmdat[-testIndexes, ]
        c.fit <- survival::survfit(survival::Surv(trainData$Y, 1 - trainData$D) ~ 1)
        cent <- testData$Y; cent[testData$D==0] <- times
        c.hat[testIndexes] <- summary(c.fit, times = cent)$surv
      }
      shudat <- data.frame(shuffle, c.hat)
      c.hat <- shudat[order(shuffle), ]$c.hat
    }else if (cen.fit == "survival.forest"){
      c.fit <- do.call(grf::survival_forest, c(list(X = cbind(X, W), Y = Y, D = 1 - D), args.nuisance))
      C.hat <- predict(c.fit, failure.times = c.fit$failure.times)$predictions
      cent <- Y; cent[D==0] <- times
      cen.times.index <- findInterval(cent, c.fit$failure.times)
      c.hat <- C.hat[cbind(1:length(Y), cen.times.index)]
    }
  }else{
    c.fit <- NULL
  }

  # create binary data
  tempdat <- data.frame(Y, D, W, m.hat, p.hat, c.hat, x.scl)
  binary.data <- tempdat[tempdat$D==1|tempdat$Y > times,]          # remove subjects who got censored before the time of interest t50
  binary.data$D[binary.data$D==1 & binary.data$Y > times] <- 0     # recode the event status for subjects who had events after t50
  binary.data <- binary.data[complete.cases(binary.data),]

  weights = 1/binary.data$c.hat     # the treatment weight is already accounted
  y.tilde = (1 - binary.data$D) - binary.data$m.hat
  x.scl = binary.data[, 7:dim(binary.data)[2]]
  foldid2 = sample(rep(seq(k.folds), length = length(binary.data$W)))

  x.scl.tilde = cbind(as.numeric(binary.data$W - binary.data$p.hat) * cbind(1, x.scl))
  x.scl.pred = cbind(1, x.scl)

  tau.fit = glmnet::cv.glmnet(as.matrix(x.scl.tilde),
                              y.tilde,
                              weights = weights,
                              foldid = foldid2,
                              alpha = 1,
                              lambda = lambda.tau,
                              penalty.factor = penalty.factor.tau, # no penalty on ATE
                              standardize = FALSE)
  tau.beta = as.vector(t(coef(tau.fit, s = lambda.choice)[-1]))
  tau.hat = as.matrix(x.scl.pred) %*% tau.beta

  ret = list(tau.fit = tau.fit,
             tau.beta = tau.beta,
             y.fit = y.fit,
             c.fit = c.fit,
             p.hat = p.hat,
             m.hat = m.hat,
             tau.hat = tau.hat)
  class(ret) <- "surv_rl_grf_lasso"
  ret
}

#' predict for surv_rl_grf_lasso
#'
#' get estimated tau(X) using the trained surv_rl_grf_lasso model
#'
#' @param object a surv_rl_grf_lasso object
#' @param newdata covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
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
#' surv.rl.grf.lasso.fit = surv_rl_grf_lasso(X, W, Y, D, times, p.hat = 0.5)
#' cate = predict(surv.rl.grf.lasso.fit)
#' cate.test = predict(surv.rl.grf.lasso.fit, X.test)
#' }
#'
#' @return A vector of predicted conditional average treatment effects
#' @export
predict.surv_rl_grf_lasso <- function(object,
                                      newdata = NULL,
                                      ...) {
  if (!is.null(newdata)) {
    newdata = sanitize_x(newdata)
    newdata.scl = scale(newdata, center = TRUE, scale = TRUE)
    newdata.scl = newdata.scl[,!is.na(colSums(newdata.scl)), drop = FALSE]
    newdata.scl.pred = cbind(1, newdata.scl)
    tau.hat = newdata.scl.pred %*% object$tau.beta
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
