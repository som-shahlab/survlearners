#' @include utils.R
#'
#' @title R-learner, implemented via glmnet (lasso)
#'
#' @description  R-learner, as proposed by Nie and Wager (2017), implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds number of folds for cross-fitting
#' @param foldid user-supplied foldid. Must have length equal to length(w). If provided, it overrides the k_folds option.
#' @param lambda_y user-supplied lambda sequence for cross validation in learning E(y|x)
#' @param lambda_w user-supplied lambda sequence for cross validation in learning E(w|x)
#' @param lambda_tau user-supplied lambda sequence for cross validation in learning the treatment effect E(y(1) - y(0) | x)
#' @param lambda_choice how to cross-validate for learning the treatment effect tau; choose from "lambda.min" or "lambda.1se"
#' @param rs whether to use the RS-learner (logical).
#' @param p_hat user-supplied estimate for E(W|X)
#' @param m_hat user-supplied estimte for E(Y|X)
#' @param penalty_factor user-supplied penalty factor, a vector of length the same as the number of covariates in x.
#' @return an rlasso object
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rlasso_fit = rlasso(x, w, y)
#' rlasso_est = predict(rlasso_fit, x)
#' }
#' @export
rlasso = function(x, w, y, D,
                  k_folds = NULL,
                  foldid = NULL,
                  lambda_y = NULL,
                  lambda_w = NULL,
                  lambda_tau = NULL,
                  lambda_choice = c("lambda.min","lambda.1se"),
                  rs = FALSE,
                  p_hat = NULL,
                  m_hat = NULL,
                  c_hat = NULL,
                  penalty_factor = NULL,
                  cf = TRUE,
                  times = NULL,
                  failure.times = NULL,
                  num.trees = 2000,
                  alpha = NULL,
                  cen_fit = "KM",
                  meta_learner = TRUE){

    input = sanitize_input(x,w,y,D)
    x = input$x
    w = input$w
    y = input$y
    D = input$D

    standardization = caret::preProcess(x, method=c("center", "scale")) # get the standardization params
    x_scl = predict(standardization, x)							 # standardize the input
    x_scl = x_scl[,!is.na(colSums(x_scl)), drop = FALSE]

    lambda_choice = match.arg(lambda_choice)

    nobs = nrow(x_scl)
    pobs = ncol(x_scl)

    if (is.null(foldid) || length(foldid) != length(w)) {

      if (!is.null(foldid) && length(foldid) != length(w)) {
        warning("supplied foldid does not have the same length ")
      }

      if (is.null(k_folds)) {
          k_folds = floor(max(3, min(10,length(w)/4)))
      }

      # fold ID for cross-validation; balance treatment assignments
      foldid = sample(rep(seq(k_folds), length = length(w)))

    }

    # penalty factor for nuisance and tau estimators
    if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
      if (!is.null(penalty_factor) && length(penalty_factor) != pobs) {
        warning("penalty_factor supplied is not of the same length as the number of columns in x after removing NA columns. Using all ones instead.")
      }
      penalty_factor_nuisance_w = rep(1, pobs)
      penalty_factor_nuisance_m = rep(1, (pobs+1))
      if (rs) {
        penalty_factor_tau = c(0, rep(1, 2 * pobs))
      }
      else {
        penalty_factor_tau = c(0, rep(1, pobs))
      }
    } else {
      penalty_factor_nuisance = penalty_factor
      if (rs) {
        penalty_factor_tau = c(0, penalty_factor, penalty_factor)
      }
      else {
        penalty_factor_tau = c(0, penalty_factor)
      }
    }

    if (is.null(p_hat)){

        if (cf){
          w_fit = glmnet::cv.glmnet(x, w,
                                    foldid = foldid,
                                    family="binomial",
                                    type.measure="deviance",
                                    keep = TRUE,
                                    lambda = lambda_w,
                                    alpha = alpha,
                                    penalty.factor = penalty_factor_nuisance_w)

          w_lambda_min = w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
          theta_hat = w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
          p_hat = 1/(1 + exp(-theta_hat))
        } else {
          w_fit = glmnet::cv.glmnet(x, w,
                                    foldid = foldid,
                                    family="binomial",
                                    type.measure="deviance",
                                    lambda = lambda_w,
                                    alpha = alpha,
                                    penalty.factor = penalty_factor_nuisance_w)

          theta_hat = as.vector(predict(w_fit, newx = x, s = w_fit$lambda.min))
          p_hat = 1/(1 + exp(-theta_hat))
        }
    }else{
      w_fit = NULL
    }

    if (is.null(m_hat)){
      if (cf){
        tempdat <- data.frame(y, D, w, x)
        foldid <- sample(rep(seq(k_folds), length = length(w)))
        survt1 <- survt0 <- rep(NA, length(w))
        for (k in 1:k_folds){
          testdat  <- tempdat[foldid==k, ]
          testdat1 <- testdat; testdat1$w <- 1
          testdat0 <- testdat; testdat0$w <- 0
          traindat <- tempdat[!foldid==k, ]
          y_fit <- glmnet::cv.glmnet(as.matrix(traindat[,3:dim(traindat)[2]]),
                                     Surv(traindat$y, traindat$D),
                                     family = "cox",
                                     nfolds = k_folds,
                                     lambda = lambda_y,
                                     alpha = alpha,
                                     penalty.factor = penalty_factor_nuisance_m)
          S0 <- base_surv(y_fit, traindat$y, traindat$D, as.matrix(traindat[,3:dim(traindat)[2]]), lambda = y_fit$lambda.min)
          survt1[foldid==k] <- pred_surv(y_fit, S0, as.matrix(cbind(testdat1[,3:dim(testdat1)[2]])), times = times, lambda = y_fit$lambda.min)
          survt0[foldid==k] <- pred_surv(y_fit, S0, as.matrix(cbind(testdat0[,3:dim(testdat0)[2]])), times = times, lambda = y_fit$lambda.min)
        }
      } else {
          traindat <- data.frame(y, D, w, x)
          testdat1 <- traindat; testdat1$w <- 1
          testdat0 <- traindat; testdat0$w <- 0
          y_fit <- glmnet::cv.glmnet(as.matrix(traindat[,3:dim(traindat)[2]]),
                                     Surv(traindat$y, traindat$D),
                                     family = "cox",
                                     nfolds = k_folds,
                                     lambda = lambda_y,
                                     alpha = alpha,
                                     penalty.factor = penalty_factor_nuisance_m)
          S0 <- base_surv(y_fit, traindat$y, traindat$D, as.matrix(traindat[,3:dim(traindat)[2]]), lambda = y_fit$lambda.min)
          survt1 <- pred_surv(y_fit, S0, as.matrix(cbind(testdat1[,3:dim(testdat1)[2]])), times = times, lambda = y_fit$lambda.min)
          survt0 <- pred_surv(y_fit, S0, as.matrix(cbind(testdat0[,3:dim(testdat0)[2]])), times = times, lambda = y_fit$lambda.min)
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

    if (is.null(failure.times)) {
      Y.grid <- sort(unique(y))
    } else {
      Y.grid <- failure.times
    }

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
      }else if (cen_fit == "lasso"){
        if (cf){
          traindat <- data.frame(y, D, w, x)
          foldid <- sample(rep(seq(k_folds), length = length(w)))
          c_fit <- glmnet::cv.glmnet(as.matrix(traindat[,3:dim(traindat)[2]]),
                                     Surv(traindat$y, 1-traindat$D),
                                     family = "cox",
                                     foldid = foldid,
                                     keep = TRUE,
                                     lambda = lambda_y,
                                     alpha = alpha,
                                     penalty.factor = penalty_factor_nuisance_m)
          c_lambda_min <- c_fit$lambda[which.min(c_fit$cvm[!is.na(colSums(c_fit$fit.preval))])]
          S0 <- base_surv(c_fit, traindat$y, 1-traindat$D, as.matrix(traindat[,3:dim(traindat)[2]]), lambda = c_lambda_min)
          cent <- rep(times, length(traindat$y))
          cent[traindat$D==1] <- traindat$y[traindat$D==1]
          c_hat <- pred_surv_preval(c_fit, S0, times = cent, lambda = c_lambda_min)
        } else {
          traindat <- data.frame(y, D, w, x)
          foldid <- sample(rep(seq(k_folds), length = length(w)))
          c_fit <- glmnet::cv.glmnet(as.matrix(traindat[,3:dim(traindat)[2]]),
                                     Surv(traindat$y, 1-traindat$D),
                                     family = "cox",
                                     foldid = foldid,
                                     lambda = lambda_y,
                                     alpha = alpha,
                                     penalty.factor = penalty_factor_nuisance_m)
          S0 <- base_surv(c_fit, traindat$y, 1-traindat$D, as.matrix(traindat[,3:dim(traindat)[2]]), lambda = c_fit$lambda.min)
          cent <- rep(times, length(traindat$y))
          cent[traindat$D==1] <- traindat$y[traindat$D==1]
          c_hat <- pred_surv(c_fit, S0, x = as.matrix(traindat[,3:dim(traindat)[2]]), times = cent, lambda = c_fit$lambda.min)
        }
      }else if (cen_fit == "survival.forest"){
        c_fit <- do.call(survival_forest, c(list(X = cbind(x, w), Y = y, D = 1 - D), args.nuisance))
        C.hat <- predict(c_fit, failure.times = Y.grid)$predictions
        cent <- y
        cent[D==0] <- times
        cen.times.index <- findInterval(cent, Y.grid)
        c_hat <- C.hat[cbind(1:length(y), cen.times.index)]
      }
    }else {
      c_fit = NULL
    }

    # use binary data
    tempdat <- data.frame(y, D, w, m_hat, p_hat, c_hat, foldid, x_scl)
    binary_data <- tempdat[tempdat$D==1|tempdat$y > times,]          # remove subjects who got censored before the time of interest t50
    binary_data$D[binary_data$D==1 & binary_data$y > times] <- 0     # recode the event status for subjects who had events after t50
    binary_data <- binary_data[complete.cases(binary_data),]

    weights = 1/binary_data$c_hat     # the treatment weight is already accounted
    y_tilde = (1 - binary_data$D) - binary_data$m_hat
    x_scl = binary_data[, 8:dim(binary_data)[2]]
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
                                  alpha = alpha,
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
    class(ret) <- "rlasso"
    ret
}


#' predict for rlasso
#'
#' get estimated tau(x) using the trained rlasso model
#'
#' @param object an rlasso object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
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
#' rlasso_fit = rlasso(x, w, y)
#' rlasso_est = predict(rlasso_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.rlasso <- function(object,
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
