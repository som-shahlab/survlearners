#' @include utils.R
#'
#' @title R-learner, implemented via xgboost (boosting)
#'
#' @description  R-learner, as proposed by Nie and Wager (2017), implemented via xgboost (boosting)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds number of folds used for cross fitting and cross validation
#' @param p_hat pre-computed estimates on E(W|X) corresponding to the input x. rboost will compute it internally if not provided.
#' @param m_hat pre-computed estimates on E(Y|X) corresponding to the input x. rboost will compute it internally if not provided.
#' @param ntrees_max the maximum number of trees to grow for xgboost
#' @param num_search_rounds the number of random sampling of hyperparameter combinations for cross validating on xgboost trees
#' @param print_every_n the number of iterations (in each iteration, a tree is grown) by which the code prints out information
#' @param early_stopping_rounds the number of rounds the test error stops decreasing by which the cross validation in finding the optimal number of trees stops
#' @param nthread the number of threads to use. The default is NULL, which uses all available threads
#' @param verbose boolean; whether to print statistic
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rboost_fit = rboost(x, w, y)
#' rboost_est = predict(rboost_fit, x)
#' }
#'
#' @export
rlasgrf = function(x, w, y, D,
                   times = NULL,
                   k_folds = NULL,
                   p_hat = NULL,
                   m_hat = NULL,
                   c_hat = NULL,
                   failure.times = NULL,
                   num.trees = 2000,
                   sample.weights = NULL,
                   clusters = NULL,
                   equalize.cluster.weights = FALSE,
                   sample.fraction = 0.5,
                   mtry = min(ceiling(sqrt(ncol(x)) + 20), ncol(x)),
                   min.node.size = 5,
                   honesty = TRUE,
                   honesty.fraction = 0.5,
                   honesty.prune.leaves = TRUE,
                   alpha = 0.05,                # splitting criteria parameter in grf
                   imbalance.penalty = 0,
                   stabilize.splits = TRUE,
                   ci.group.size = 2,
                   tune.parameters = "none",
                   compute.oob.predictions = TRUE,
                   num.threads = NULL,
                   seed = runif(1, 0, .Machine$integer.max),
                   lambda_tau = NULL,           # for lasso
                   rs = FALSE,
                   penalty_factor = NULL,
                   lambda_choice = c("lambda.min","lambda.1se"),
                   cen_fit = "KM",
                   meta_learner = TRUE,
                   verbose = FALSE){


  input = sanitize_input(x,w,y,D)
  x = input$x
  w = as.numeric(input$w)
  y = input$y
  D = input$D
  nobs = nrow(x)
  pobs = ncol(x)

  standardization = caret::preProcess(x, method=c("center", "scale"))  # get the standardization params
  x_scl = predict(standardization, x)				                    			 # standardize the input
  x_scl = x_scl[,!is.na(colSums(x_scl)), drop = FALSE]

  lambda_choice = match.arg(lambda_choice)

  # penalty factor for tau estimator
  if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
    if (!is.null(penalty_factor) && length(penalty_factor) != pobs) {
      warning("penalty_factor supplied is not of the same length as the number of columns in x after removing NA columns. Using all ones instead.")
    }
    if (rs) {
      penalty_factor_tau = c(0, rep(1, 2 * pobs))
    }
    else {
      penalty_factor_tau = c(0, rep(1, pobs))
    }
  } else {
    if (rs) {
      penalty_factor_tau = c(0, penalty_factor, penalty_factor)
    }
    else {
      penalty_factor_tau = c(0, penalty_factor)
    }
  }

  if (is.null(p_hat)){
    w_fit <- regression_forest(x, w, num.trees = max(50, num.trees / 4),
                               sample.weights = sample.weights, clusters = clusters,
                               equalize.cluster.weights = equalize.cluster.weights,
                               sample.fraction = sample.fraction, mtry = mtry,
                               min.node.size = 5, honesty = TRUE,
                               honesty.fraction = 0.5, honesty.prune.leaves = TRUE,
                               alpha = alpha, imbalance.penalty = imbalance.penalty,
                               ci.group.size = 1, compute.oob.predictions = TRUE,
                               num.threads = num.threads, seed = seed)
    p_hat <- predict(w_fit)$predictions
  }else if (length(p_hat) == 1) {
    w_fit = NULL
    p_hat <- rep(p_hat, nrow(x))
  }else if (length(p_hat) != nrow(x)){
    stop("p_hat has incorrect length.")
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

  if (is.null(m_hat)){
    y_fit <- do.call(survival_forest, c(list(X = cbind(x, w), Y = y, D = D), args.nuisance))
    y_fit[["X.orig"]][, ncol(x) + 1] <- rep(1, nrow(x))
    S1.hat <- predict(y_fit)$predictions
    y_fit[["X.orig"]][, ncol(x) + 1] <- rep(0, nrow(x))
    S0.hat <- predict(y_fit)$predictions
    y_fit[["X.orig"]][, ncol(x) + 1] <- w

    times.index <- findInterval(times, y_fit$failure.times)
    surf1 <- S1.hat[, times.index]
    surf0 <- S0.hat[, times.index]
    m_hat  <- p_hat * surf1 + (1 - p_hat) * surf0
  }else {
    y_fit = NULL
  }

  if (is.null(failure.times)) {
    Y.grid <- sort(unique(y))
  } else {
    Y.grid <- failure.times
  }

  args.nuisance$compute.oob.predictions <- TRUE
  if (is.null(c_hat)){
    if(cen_fit == "KM"){
      traindat <- data.frame(Y = y, D = D)
      shuffle <- sample(nrow(traindat))
      kmdat <- traindat[shuffle,]
      folds <- cut(seq(1, nrow(kmdat)), breaks=10, labels=FALSE)
      c_hat <- rep(NA, nrow(kmdat))
      for(z in 1:10){
        testIndexes <- which(folds==z, arr.ind=TRUE)
        testData <- kmdat[testIndexes, ]
        trainData <- kmdat[-testIndexes, ]
        c_fit <- survfit(Surv(trainData$Y, 1 - trainData$D) ~ 1)
        cent <- testData$Y
        cent[testData$D==0] <- times
        c_hat[testIndexes] <- summary(c_fit, times = cent)$surv
      }
      shudat <- data.frame(shuffle, c_hat)
      c_hat <- shudat[order(shuffle), ]$c_hat
    }else if (cen_fit == "survival.forest"){
      c_fit <- do.call(survival_forest, c(list(X = cbind(x, w), Y = y, D = 1 - D), args.nuisance))
      C.hat <- predict(c_fit, failure.times = Y.grid)$predictions
      cent <- y
      cent[D==0] <- times
      cen.times.index <- findInterval(cent, Y.grid)
      c_hat <- C.hat[cbind(1:length(y), cen.times.index)]
    }
  }else{
    c_fit <- NULL
  }

  # create binary data
  tempdat <- data.frame(y, D, w, m_hat, p_hat, c_hat, x_scl)
  binary_data <- tempdat[tempdat$D==1|tempdat$y > times,]          # remove subjects who got censored before the time of interest t50
  binary_data$D[binary_data$D==1 & binary_data$y > times] <- 0     # recode the event status for subjects who had events after t50
  binary_data <- binary_data[complete.cases(binary_data),]

  weights = 1/binary_data$c_hat     # the treatment weight is already accounted
  y_tilde = (1 - binary_data$D) - binary_data$m_hat
  x_scl = binary_data[, 7:dim(binary_data)[2]]
  foldid2 = sample(rep(seq(k_folds), length = length(binary_data$w)))

  if (rs){
    x_scl_tilde = cbind(as.numeric(binary_data$w - binary_data$p_hat) * cbind(1, x_scl), x_scl)
    x_scl_pred = cbind(1, x_scl, x_scl * 0)
  }else{
    x_scl_tilde = cbind(as.numeric(binary_data$w - binary_data$p_hat) * cbind(1, x_scl))
    x_scl_pred = cbind(1, x_scl)
  }

  if (meta_learner){
    tau_fit = glmnet::cv.glmnet(as.matrix(x_scl_tilde),
                                y_tilde,
                                weights = weights,
                                foldid = foldid2,
                                alpha = 1,
                                lambda = lambda_tau,
                                penalty.factor = penalty_factor_tau, # no penalty on ATE
                                standardize = FALSE)
    tau_beta = as.vector(t(coef(tau_fit, s = lambda_choice)[-1]))
    tau_hat = as.matrix(x_scl_pred) %*% tau_beta
  }else{
    dat = data.frame(y_tilde, x_scl_tilde)
    tau_fit = glm(y_tilde ~ .,
                  family = "gaussian",
                  weights = weights,
                  data = dat)
    tau_beta = as.vector(t(tau_fit$coefficients[-1]))
    tau_hat = as.matrix(x_scl_pred) %*% tau_beta
  }

  ret = list(tau_fit = tau_fit,
             tau_beta = tau_beta,
             w_fit = w_fit,
             y_fit = y_fit,
             c_fit = c_fit,
             p_hat = p_hat,
             m_hat = m_hat,
             tau_hat = tau_hat,
             rs = rs,
             standardization = standardization)
  class(ret) <- "rlasgrf"
  ret
}

#' predict for rgrf
#'
#' get estimated tau(x) using the trained rgrf model
#'
#' @param object a rgrf object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param tau_only if set to TRUE, onlly return prediction on tau. Otherwise, return a list including prediction on tau, propensity score, and baseline main effect.
#' @param ... additional arguments (currently not used)
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rboost_fit = rboost(x, w, y)
#' rboost_est = predict(rboost_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.rlasgrf <- function(object,
                            newx = NULL,
                            ...) {
  if (!is.null(newx)) {

    newx = sanitize_x(newx)
    newx_scl = predict(object$standardization, newx) # standardize the new data using the same standardization as with the training data
    newx_scl = newx_scl[,!is.na(colSums(newx_scl)), drop = FALSE]

    if (object$rs){
      newx_scl_pred = cbind(1, newx_scl, newx_scl * 0)
    }
    else{
      newx_scl_pred = cbind(1, newx_scl)
    }
    tau_hat = newx_scl_pred %*% object$tau_beta
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
