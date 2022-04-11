#' @title R-learner of lasso
#'
#' @description  R-learner, implemented via glmnet (lasso)
#'
#' @param X The baseline covariates
#' @param W The treatment variable (0 or 1)
#' @param Y The follow-up time
#' @param D The event indicator
#' @param k_folds Number of folds for cross validation
#' @param foldid User-supplied foldid. Must have length equal to length(W). If provided, it overrides the k_folds option.
#' @param lambda_y User-supplied lambda sequence for cross validation in the outcome model
#' @param lambda_tau User-supplied lambda sequence for cross validation in the cate model
#' @param lambda_choice How to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param p_hat Propensity score
#' @param m_hat Conditional mean outcome E(Y|X)
#' @param c_hat Censoring weights
#' @param penalty_factor User-supplied penalty factor, must be of length the same as number of features in X
#' @param times The prediction time of interest
#' @param failure.times A vector of event times to fit the survival curve at.
#' @param num.trees Number of trees grown in the forest
#' @param alpha Imbalance tuning parameter for a split in a tree
#' @param cen_fit The choice of model fitting for censoring
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
#' X.test <- matrix(rnorm(n * p), n, p)
#'
#' surv_rl_lasso_fit = surv_rl_lasso(X, W, Y, D, times, p_hat = 0.5)
#' cate = predict(surv_rl_lasso_fit)
#' cate.test = predict(surv_rl_lasso_fit, X.test)
#' }
#' @return a surv_rl_lasso object
#' @export
surv_rl_lasso = function(X, W, Y, D,
                         times = NULL,
                         k_folds = 10,
                         foldid = NULL,
                         lambda_y = NULL,
                         lambda_tau = NULL,
                         lambda_choice = "lambda.min",
                         p_hat = NULL,
                         m_hat = NULL,
                         c_hat = NULL,
                         penalty_factor = NULL,
                         failure.times = NULL,
                         num.trees = 2000,
                         alpha = 0.05,
                         cen_fit = "KM"){

    input = sanitize_input(X, W, Y, D)
    X = input$X
    W = input$W
    Y = input$Y
    D = input$D

    x_scl = scale(X, center = TRUE, scale = TRUE)
    x_scl = x_scl[,!is.na(colSums(x_scl)), drop = FALSE]

    nobs = nrow(x_scl)
    pobs = ncol(x_scl)

    if (is.null(foldid) || length(foldid) != length(W)) {

      if (!is.null(foldid) && length(foldid) != length(W)) {
        warning("supplied foldid does not have the same length ")
      }

      if (is.null(k_folds)) {
          k_folds = floor(max(3, min(10,length(W)/4)))
      }

      # fold ID for cross-validation; balance treatment assignments
      foldid = sample(rep(seq(k_folds), length = length(W)))

    }

    # penalty factor for nuisance and tau estimators
    if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
      if (!is.null(penalty_factor) && length(penalty_factor) != pobs) {
        warning("penalty_factor supplied is not of the same length as the number of columns in X after removing NA columns. Using all ones instead.")
      }
      penalty_factor_nuisance_w = rep(1, pobs)
      penalty_factor_nuisance_m = rep(1, (pobs+1))
      penalty_factor_tau = c(0, rep(1, pobs))
    } else {
      penalty_factor_nuisance = penalty_factor
      penalty_factor_tau = c(0, penalty_factor)
    }

    if (is.null(p_hat)){
      stop("propensity score needs to be supplied")
    }else if (length(p_hat) == 1) {
      p_hat <- rep(p_hat, nrow(X))
    }else if (length(p_hat) != nrow(X)){
      stop("p_hat has incorrect length.")
    }

    if (is.null(m_hat)){
    foldid <- sample(rep(seq(k_folds), length = length(W)))
    survt1 <- survt0 <- rep(NA, length(W))
    for (k in 1:k_folds){
      XW <- as.matrix(data.frame(W[!foldid==k], X[!foldid==k, ]))
      y_fit <- glmnet::cv.glmnet(XW,
                                 survival::Surv(Y[!foldid==k], D[!foldid==k]),
                                 family = "cox",
                                 nfolds = 10,
                                 lambda = lambda_y,
                                 alpha = 1,
                                 penalty.factor = penalty_factor_nuisance_m)
      S0 <- base_surv(y_fit, Y[!foldid==k], D[!foldid==k], XW, lambda = y_fit$lambda.min)
      survt1[foldid==k] <- pred_surv(y_fit, S0, cbind(rep(1, length(W[foldid==k])), X[foldid==k, ]), times = times, lambda = y_fit$lambda.min)
      survt0[foldid==k] <- pred_surv(y_fit, S0, cbind(rep(0, length(W[foldid==k])), X[foldid==k, ]), times = times, lambda = y_fit$lambda.min)
    }
    m_hat  <- p_hat * survt1 + (1 - p_hat) * survt0
    }else {
      y_fit = NULL
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

    if (is.null(c_hat)){
    if(cen_fit == "KM"){
      traindat <- data.frame(Y = Y, D = D)
      shuffle <- sample(nrow(traindat))
      kmdat <- traindat[shuffle,]
      folds <- cut(seq(1, nrow(kmdat)), breaks=10, labels=FALSE)
      c_hat <- rep(NA, nrow(kmdat))
      for(z in 1:10){
        testIndexes <- which(folds==z, arr.ind=TRUE)
        testData <- kmdat[testIndexes, ]
        trainData <- kmdat[-testIndexes, ]
        c_fit <- survival::survfit(survival::Surv(trainData$Y, 1 - trainData$D) ~ 1)
        cent <- testData$Y; cent[testData$D==0] <- times
        c_hat[testIndexes] <- summary(c_fit, times = cent)$surv
      }
      shudat <- data.frame(shuffle, c_hat)
      c_hat <- shudat[order(shuffle), ]$c_hat
    }else if (cen_fit == "survival.forest"){
      cc_fit <- do.call(grf::survival_forest, c(list(X = cbind(X, W), Y = Y, D = 1 - D), args.nuisance))
      C.hat <- predict(c_fit, failure.times = c_fit$failure.times)$predictions
      cent <- Y; cent[D==0] <- times
      cen.times.index <- findInterval(cent, c_fit$failure.times)
      c_hat <- C.hat[cbind(1:length(Y), cen.times.index)]
     }
    }else {
      c_fit = NULL
    }

    # use binary data
    tempdat <- data.frame(Y, D, W, m_hat, p_hat, c_hat, foldid, x_scl)
    binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          # remove subjects who got censored before the time of interest t50
    binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     # recode the event status for subjects who had events after t50
    binary_data <- binary_data[complete.cases(binary_data),]

    weights = 1/binary_data$c_hat
    y_tilde = (1 - binary_data$D) - binary_data$m_hat
    x_scl = binary_data[, 8:dim(binary_data)[2]]
    foldid2 = sample(rep(seq(k_folds), length = length(binary_data$W)))

    x_scl_tilde = cbind(as.numeric(binary_data$W - binary_data$p_hat) * cbind(1, x_scl))
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
               c_hat = c_hat,
               tau_hat = tau_hat)
    class(ret) <- "surv_rl_lasso"
    ret
}


#' predict for surv_rl_lasso
#'
#' get estimated tau(X) using the trained surv_rl_lasso model
#'
#' @param object An surv_rl_lasso object
#' @param newx Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
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
#' X.test <- matrix(rnorm(n * p), n, p)
#'
#' surv_rl_lasso_fit = surv_rl_lasso(X, W, Y, D, times, p_hat = 0.5)
#' cate = predict(surv_rl_lasso_fit)
#' cate.test = predict(surv_rl_lasso_fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_rl_lasso <- function(object,
                                  newx = NULL,
                                  ...) {
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    newx_scl = scale(newx, center = TRUE, scale = TRUE)
    newx_scl = newx_scl[,!is.na(colSums(newx_scl)), drop = FALSE]
    newx_scl_pred = cbind(1, newx_scl)
    tau_hat = newx_scl_pred %*% object$tau_beta
  }else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
