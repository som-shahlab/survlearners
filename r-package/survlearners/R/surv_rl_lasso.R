#' @title R-learner of lasso
#'
#' @description  R-learner, implemented via glmnet (lasso)
#'
#' @param X The baseline covariates
#' @param W The treatment variable (0 or 1)
#' @param Y The follow-up time
#' @param D The event indicator
#' @param k.folds Number of folds for cross validation
#' @param foldid User-supplied foldid. Must have length equal to length(W). If provided, it overrides the k.folds option.
#' @param lambda.y User-supplied lambda sequence for cross validation in the outcome model
#' @param lambda.tau User-supplied lambda sequence for cross validation in the cate model
#' @param lambda.choice How to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param p.hat Propensity score
#' @param m.hat Conditional mean outcome E(Y|X)
#' @param c.hat Censoring weights
#' @param penalty.factor User-supplied penalty factor, must be of length the same as number of features in X
#' @param times The prediction time of interest
#' @param failure.times A vector of event times to fit the survival curve at.
#' @param num.trees Number of trees grown in the forest
#' @param alpha Imbalance tuning parameter for a split in a tree
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
#' surv.rl.lasso.fit = surv_rl_lasso(X, W, Y, D, times, p.hat = 0.5)
#' cate = predict(surv.rl.lasso.fit)
#' cate.test = predict(surv.rl.lasso.fit, X.test)
#' }
#' @return a surv_rl_lasso object
#' @export
surv_rl_lasso = function(X, W, Y, D,
                         times = NULL,
                         k.folds = 10,
                         foldid = NULL,
                         lambda.y = NULL,
                         lambda.tau = NULL,
                         lambda.choice = "lambda.min",
                         p.hat = NULL,
                         m.hat = NULL,
                         c.hat = NULL,
                         penalty.factor = NULL,
                         failure.times = NULL,
                         num.trees = 2000,
                         alpha = 0.05,
                         cen.fit = "KM"){

    input = sanitize_input(X, W, Y, D)
    X = input$X
    W = input$W
    Y = input$Y
    D = input$D

    x.scl = scale(X, center = TRUE, scale = TRUE)
    x.scl = x.scl[,!is.na(colSums(x.scl)), drop = FALSE]

    nobs = nrow(x.scl)
    pobs = ncol(x.scl)

    if (is.null(foldid) || length(foldid) != length(W)) {

      if (!is.null(foldid) && length(foldid) != length(W)) {
        warning("supplied foldid does not have the same length ")
      }

      if (is.null(k.folds)) {
          k.folds = floor(max(3, min(10,length(W)/4)))
      }

      # fold ID for cross-validation; balance treatment assignments
      foldid = sample(rep(seq(k.folds), length = length(W)))

    }

    # penalty factor for nuisance and tau estimators
    if (is.null(penalty.factor) || (length(penalty.factor) != pobs)) {
      if (!is.null(penalty.factor) && length(penalty.factor) != pobs) {
        warning("penalty.factor supplied is not of the same length as the number of columns in X after removing NA columns. Using all ones instead.")
      }
      penalty.factor.nuisance.w = rep(1, pobs)
      penalty.factor.nuisance.m = rep(1, (pobs+1))
      penalty.factor.tau = c(0, rep(1, pobs))
    } else {
      penalty.factor.nuisance = penalty.factor
      penalty.factor.tau = c(0, penalty.factor)
    }

    if (is.null(p.hat)){
      stop("propensity score needs to be supplied")
    }else if (length(p.hat) == 1) {
      p.hat <- rep(p.hat, nrow(X))
    }else if (length(p.hat) != nrow(X)){
      stop("p.hat has incorrect length.")
    }

    if (is.null(m.hat)){
    foldid <- sample(rep(seq(k.folds), length = length(W)))
    survt1 <- survt0 <- rep(NA, length(W))
    for (k in 1:k.folds){
      XW <- as.matrix(data.frame(W[!foldid==k], X[!foldid==k, ]))
      y.fit <- glmnet::cv.glmnet(XW,
                                 survival::Surv(Y[!foldid==k], D[!foldid==k]),
                                 family = "cox",
                                 nfolds = 10,
                                 lambda = lambda.y,
                                 alpha = 1,
                                 penalty.factor = penalty.factor.nuisance.m)
      S0 <- base_surv(y.fit, Y[!foldid==k], D[!foldid==k], XW, lambda = y.fit$lambda.min)
      survt1[foldid==k] <- pred_surv(y.fit, S0, cbind(rep(1, length(W[foldid==k])), X[foldid==k, ]), times = times, lambda = y.fit$lambda.min)
      survt0[foldid==k] <- pred_surv(y.fit, S0, cbind(rep(0, length(W[foldid==k])), X[foldid==k, ]), times = times, lambda = y.fit$lambda.min)
    }
    m.hat  <- p.hat * survt1 + (1 - p.hat) * survt0
    }else {
      y.fit = NULL
    }

    args.nuisance <- list(failure.times = failure.times,
                          num.trees = max(50, num.trees / 4),
                          min.node.size = 15,
                          honesty = TRUE,
                          honesty.fraction = 0.5,
                          honesty.prune.leaves = TRUE,
                          alpha = alpha,
                          prediction.type = "Nelson-Aalen",
                          compute.oob.predictions = TRUE)

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
    }else {
      c.fit = NULL
    }

    # use binary data
    tempdat <- data.frame(Y, D, W, m.hat, p.hat, c.hat, foldid, x.scl)
    binary.data <- tempdat[tempdat$D==1|tempdat$Y > times,]          # remove subjects who got censored before the time of interest t50
    binary.data$D[binary.data$D==1 & binary.data$Y > times] <- 0     # recode the event status for subjects who had events after t50
    binary.data <- binary.data[complete.cases(binary.data),]

    weights = 1/binary.data$c.hat
    y.tilde = (1 - binary.data$D) - binary.data$m.hat
    x.scl = binary.data[, 8:dim(binary.data)[2]]
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
               c.hat = c.hat,
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
#' surv.rl.lasso.fit = surv_rl_lasso(X, W, Y, D, times, p.hat = 0.5)
#' cate = predict(surv.rl.lasso.fit)
#' cate.test = predict(surv.rl.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_rl_lasso <- function(object,
                                  newdata = NULL,
                                  ...) {
  if (!is.null(newdata)) {
    newdata = sanitize_x(newdata)
    newdata.scl = scale(newdata, center = TRUE, scale = TRUE)
    newdata.scl = newdata.scl[,!is.na(colSums(newdata.scl)), drop = FALSE]
    newdata.scl.pred = cbind(1, newdata.scl)
    tau.hat = newdata.scl.pred %*% object$tau.beta
  }else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
