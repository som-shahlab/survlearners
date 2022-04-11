#' @title T-learner of lasso
#'
#' @description  T-learner, implemented via glmnet (lasso)
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
#' @param newX The test data set (covariates only)
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
#'
#' surv_tl_lasso_fit = surv_tl_lasso(X, W, Y, D, times)
#' cate = predict(surv_tl_lasso_fit)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_tl_lasso <- function(X, W, Y, D, times, newX = NULL){

  # Model for W = 1
  foldid <- sample(rep(seq(10), length = length(Y[W==1])))
  x1 <- as.matrix(data.frame(X[W==1, ]))
  lasso_fit1 <- glmnet::cv.glmnet(x1,
                                  survival::Surv(Y[W==1], D[W==1]),
                                  family = "cox",
                                  alpha = 1,
                                  foldid = foldid)

  bsurv1 <- base_surv(fit = lasso_fit1,
                    Y = Y[W==1],
                    D = D[W==1],
                    X = x1,
                    lambda = lasso_fit1$lambda.min)

  surf1 <- pred_surv(fit = lasso_fit1,
                      S0 = bsurv1,
                       X = X,
                       times = times,
                       lambda = lasso_fit1$lambda.min)


  # Model for W = 0
  foldid <- sample(rep(seq(10), length = length(Y[W==0])))
  x0 <- as.matrix(data.frame(X[W==0, ]))
  lasso_fit0 <- glmnet::cv.glmnet(x0,
                                  survival::Surv(Y[W==0], D[W==0]),
                                  family = "cox",
                                  alpha = 1,
                                  foldid = foldid)

  bsurv0 <- base_surv(fit = lasso_fit0,
                    Y = Y[W==0],
                    D = D[W==0],
                    X = x0,
                    lambda = lasso_fit0$lambda.min)

  surf0 <- pred_surv(fit = lasso_fit0,
                       S0 = bsurv0,
                       X = X,
                       times = times,
                       lambda = lasso_fit0$lambda.min)

  pred_T_lasso <- surf1 - surf0

  ret <- list(fit1 = lasso_fit1,
              fit0 = lasso_fit0,
              bsurv1 = bsurv1,
              bsurv0 = bsurv0,
              tau = pred_T_lasso,
              times = times)
  class(ret) <- 'surv_tl_lasso'
  ret
}

#' predict for surv_tl_lasso
#'
#' get estimated tau(X) using the trained surv_tl_lasso model
#'
#' @param object An surv_tl_lasso object
#' @param newx Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param times The prediction time of interest
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
#'
#' surv_tl_lasso_fit = surv_tl_lasso(X, W, Y, D, times)
#' cate = predict(surv_tl_lasso_fit)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_tl_lasso = function(object,
                                 newx = NULL,
                                 times = NULL,
                                 ...) {
  if(is.null(newx)){
    return(object$tau)
  }else{
    if(is.null(times)){
      surf1 <- pred_surv(fit = object$fit1,
                         S0 = object$bsurv1,
                         X = newX,
                         times = object$times,
                         lambda = object$fit1$lambda.min)

      surf0 <- pred_surv(fit = object$fit0,
                         S0 = object$bsurv0,
                         X = newX,
                         times = object$times,
                         lambda = object$fit0$lambda.min)
    }else{
      surf1 <- pred_surv(fit = object$fit1,
                         S0 = object$bsurv1,
                         X = newX,
                         times = times,
                         lambda = object$fit1$lambda.min)

      surf0 <- pred_surv(fit = object$fit0,
                         S0 = object$bsurv0,
                         X = newX,
                         times = times,
                         lambda = object$fit0$lambda.min)
    }
    return(surf1 - surf0)
  }
}
