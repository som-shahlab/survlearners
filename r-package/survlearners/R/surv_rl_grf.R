#' @title R-learner of grf
#'
#' @description  R-learner, implemented via the grf package
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
#' surv.rl.grf.fit = surv_rl_grf(X, W, Y, D, times, p.hat = 0.5)
#' cate = predict(surv.rl.grf.fit)
#' cate.test = predict(surv.rl.grf.fit, X.test)
#' }
#' @return a surv_rl_grf_fit object
#' @export
surv_rl_grf = function(X, W, Y, D,
                       times = NULL,
                       k.folds = NULL,
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
                       cen.fit = "KM",
                       verbose = FALSE){

  input = sanitize_input(X,W,Y,D)
  X = input$X
  W = as.numeric(input$W)
  Y = input$Y
  D = input$D
  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(k.folds)) {
    k.folds = floor(max(3, min(10,length(Y)/4)))
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
  tempdata <- data.frame(Y, D, W, m.hat, p.hat, c.hat, X)
  binary.data <- tempdata[tempdata$D==1|tempdata$Y > times,]       # remove subjects who got censored before the time of interest t50
  binary.data$D[binary.data$D==1 & binary.data$Y > times] <- 0     # recode the event status for subjects who had events after t50
  binary.data <- binary.data[complete.cases(binary.data),]

  y.tilde <- (1 - binary.data$D) - binary.data$m.hat
  w.tilde <-  binary.data$W - binary.data$p.hat
  pseudo.outcome <- y.tilde/w.tilde
  weights <- w.tilde^2/binary.data$c.hat

  tau.fit <- grf::regression_forest(binary.data[,7:dim(binary.data)[2]],
                                    pseudo.outcome,
                                    sample.weights = weights,
                                    num.trees = num.trees,
                                    clusters = clusters,
                                    sample.fraction = sample.fraction,
                                    mtry = mtry,
                                    min.node.size = min.node.size,
                                    honesty = honesty,
                                    honesty.fraction = honesty.fraction,
                                    honesty.prune.leaves = honesty.prune.leaves,
                                    imbalance.penalty = imbalance.penalty,
                                    ci.group.size = ci.group.size,
                                    compute.oob.predictions = compute.oob.predictions,
                                    num.threads = num.threads,
                                    seed = seed)

  ret <- list(tau.fit = tau.fit,
              pseudo.outcome = pseudo.outcome,
              weights = weights,
              y.fit = y.fit,
              c.fit = c.fit,
              p.hat = p.hat,
              m.hat = m.hat,
              c.hat = c.hat)
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
#' surv.rl.grf.fit = surv_rl_grf(X, W, Y, D, times, p.hat = 0.5)
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
