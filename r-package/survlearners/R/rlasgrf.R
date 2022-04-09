#' @title R-learner of grf and lasso
#'
#' @description  R-learner, implemented via the grf package for nuisance parameter estimation and lasso for target parameter
#'
#' @param x The baseline covariates
#' @param w The treatment variable (0 or 1)
#' @param y The follow-up time
#' @param D The event indicator
#' @param k_folds Number of folds for cross validation
#' @param p_hat Propensity score
#' @param m_hat Conditional mean outcome E[Y|X]
#' @param c_hat Censoring weights
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
#' @param lambda_tau User-supplied lambda sequence for cross validation in the cate model
#' @param lambda_choice How to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty_factor User-supplied penalty factor, must be of length the same as number of features in x
#' @param cen_fit The choice of model fitting for censoring
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' times <- 0.2
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
#'
#' rlasgrf_fit = rlasgrf(X, W, Y, D, times)
#' rlasgrf_cate = predict(rlasgrf_fit, X)
#' }
#' @return a rlasgrf object
#' @export
rlasgrf = function(x, w, y, D,
                   times = NULL,
                   k_folds = 10,
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
                   alpha = 0.05,
                   imbalance.penalty = 0,
                   stabilize.splits = TRUE,
                   ci.group.size = 2,
                   tune.parameters = "none",
                   compute.oob.predictions = TRUE,
                   num.threads = NULL,
                   seed = runif(1, 0, .Machine$integer.max),
                   lambda_tau = NULL,
                   lambda_choice = "lambda.min",
                   penalty_factor = NULL,
                   cen_fit = "KM",
                   verbose = FALSE){

  input = sanitize_input(x,w,y,D)
  x = input$x
  w = as.numeric(input$w)
  y = input$y
  D = input$D
  nobs = nrow(x)
  pobs = ncol(x)

  x_scl = scale(x, center = TRUE, scale = TRUE)
  x_scl = x_scl[,!is.na(colSums(x_scl)), drop = FALSE]

  # penalty factor for tau estimator
  if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
    if (!is.null(penalty_factor) && length(penalty_factor) != pobs) {
      warning("penalty_factor supplied is not of the same length as the number of columns in x after removing NA columns. Using all ones instead.")
    }
    penalty_factor_tau = c(0, rep(1, pobs))
  } else {
    penalty_factor_tau = c(0, penalty_factor)
  }

  if (is.null(p_hat)){
    stop("propensity score needs to be supplied")
  }else if (length(p_hat) == 1) {
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
        cent <- testData$Y; cent[testData$D==0] <- times
        c_hat[testIndexes] <- summary(c_fit, times = cent)$surv
      }
      shudat <- data.frame(shuffle, c_hat)
      c_hat <- shudat[order(shuffle), ]$c_hat
    }else if (cen_fit == "survival.forest"){
      c_fit <- do.call(survival_forest, c(list(X = cbind(x, w), Y = y, D = 1 - D), args.nuisance))
      C.hat <- predict(c_fit, failure.times = c_fit$failure.times)$predictions
      cent <- y; cent[D==0] <- times
      cen.times.index <- findInterval(cent, c_fit$failure.times)
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

  x_scl_tilde = cbind(as.numeric(binary_data$w - binary_data$p_hat) * cbind(1, x_scl))
  x_scl_pred = cbind(1, x_scl)

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

  ret = list(tau_fit = tau_fit,
             tau_beta = tau_beta,
             y_fit = y_fit,
             c_fit = c_fit,
             p_hat = p_hat,
             m_hat = m_hat,
             tau_hat = tau_hat)
  class(ret) <- "rlasgrf"
  ret
}

#' predict for rlasgrf
#'
#' get estimated tau(x) using the trained rlasgrf model
#'
#' @param object a rlasgrf object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param ... additional arguments (currently not used)
#'
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' times <- 0.2
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
#'
#' rlasgrf_fit = rlasgrf(X, W, Y, D, times)
#' rlasgrf_cate = predict(rlasgrf_fit, X)
#' }
#'
#' @return A vector of predicted conditional average treatment effects
#' @export
predict.rlasgrf <- function(object,
                            newx = NULL,
                            ...) {
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    newx_scl = scale(newx, center = TRUE, scale = TRUE)
    newx_scl = newx_scl[,!is.na(colSums(newx_scl)), drop = FALSE]
    newx_scl_pred = cbind(1, newx_scl)
    tau_hat = newx_scl_pred %*% object$tau_beta
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
