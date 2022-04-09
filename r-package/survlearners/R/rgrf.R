#' @title R-learner of grf
#'
#' @description  R-learner, implemented via the grf package
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
#' rgrf_fit = rgrf(X, W, Y, D, times)
#' rgrf_cate = predict(rgrf_fit, X)
#' }
#' @return a rgrf object
#' @export
rgrf = function(x, w, y, D,
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
                alpha = 0.05,
                imbalance.penalty = 0,
                stabilize.splits = TRUE,
                ci.group.size = 2,
                tune.parameters = "none",
                compute.oob.predictions = TRUE,
                num.threads = NULL,
                seed = runif(1, 0, .Machine$integer.max),
                cen_fit = "KM",
                verbose = FALSE){

  input = sanitize_input(x,w,y,D)
  x = input$x
  w = as.numeric(input$w)
  y = input$y
  D = input$D
  nobs = nrow(x)
  pobs = ncol(x)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(y)/4)))
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
  tempdata <- data.frame(y, D, w, m_hat, p_hat, c_hat, x)
  binary_data <- tempdata[tempdata$D==1|tempdata$y > times,]       # remove subjects who got censored before the time of interest t50
  binary_data$D[binary_data$D==1 & binary_data$y > times] <- 0     # recode the event status for subjects who had events after t50
  binary_data <- binary_data[complete.cases(binary_data),]

  y_tilde <- (1 - binary_data$D) - binary_data$m_hat
  w_tilde <-  binary_data$w - binary_data$p_hat
  pseudo_outcome <- y_tilde/w_tilde
  weights <- w_tilde^2/binary_data$c_hat

  tau_fit <- grf::regression_forest(binary_data[,7:dim(binary_data)[2]],
                                    pseudo_outcome,
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

  ret <- list(tau_fit = tau_fit,
              pseudo_outcome = pseudo_outcome,
              weights = weights,
              y_fit = y_fit,
              c_fit = c_fit,
              p_hat = p_hat,
              m_hat = m_hat,
              c_hat = c_hat)
  class(ret) <- "rgrf"
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
#' rgrf_fit = rgrf(X, W, Y, D, times)
#' rgrf_cate = predict(rgrf_fit, X)
#' }
#'
#' @return A vector of predicted conditional average treatment effects
#' @export
predict.rgrf <- function(object,
                         newx = NULL,
                         tau_only = TRUE,
                          ...) {
  if (!is.null(newx)){
    newx = sanitize_x(newx)
  }
  if (tau_only) {
    return(predict(object$tau_fit, newx)$predictions)
  } else {
    tau <- predict(object$tau_fit, newx)$predictions
    e = predict(object$w_fit, newx)$predictions
    m = predict(object$y_fit, newx)$predictions
    mu1 = m + (1-e) * tau
    mu0 = m - e * tau
    return(list(tau=tau, e=e, m=m, mu1 = mu1, mu0 = mu0))
  }
}
