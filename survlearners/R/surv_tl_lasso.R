#' @title T-learner of lasso
#'
#' @description  T-learner, implemented via glmnet (lasso)
#'
#' @param data The training data set
#' @param data.test The testing data set
#' @param times The prediction time of interest
#' @examples
#' \dontrun{
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
#' data <- list(X = X, W = W, Y = Y, D = D)
#' data.test <- list(X = X, W = W, Y = Y, D = D)
#'
#' cate = surv_tl_lasso(data, data.test, times)
#' }
#' @return A vector of estimated conditional average treatment effects
#' @export
surv_tl_lasso <- function(data, data.test, times){

  traindat <- data.frame(Y = data$Y, D = data$D, W = data$W, data$X)
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]

  # Model for W = 1
  foldid <- sample(rep(seq(10), length = length(traindat1$Y)))
  lasso_fit1 <- glmnet::cv.glmnet(as.matrix(traindat1[,3:dim(traindat1)[2]]),
                                  Surv(traindat1$Y, traindat1$D),
                                  family = "cox",
                                  alpha = 1,
                                  foldid = foldid)

  S0 <- base_surv(fit = lasso_fit1,
                    Y = traindat1$Y,
                    D = traindat1$D,
                    x = as.matrix(traindat1[,3:dim(traindat1)[2]]),
                    lambda = lasso_fit1$lambda.min)

  surf1 <- pred_surv(fit = lasso_fit1,
                      S0 = S0,
                       x = data.test$X,
                       times = times,
                       lambda = lasso_fit1$lambda.min)


  # Model for W = 0
  foldid <- sample(rep(seq(10), length = length(traindat0$Y)))
  lasso_fit0 <- glmnet::cv.glmnet(as.matrix(traindat0[,3:dim(traindat0)[2]]),
                          Surv(traindat0$Y, traindat0$D),
                          family = "cox",
                          alpha = 1,
                          foldid = foldid)

  S0 <- base_surv(fit = lasso_fit0,
                    Y = traindat0$Y,
                    D = traindat0$D,
                    x = as.matrix(traindat0[,3:dim(traindat0)[2]]),
                    lambda = lasso_fit0$lambda.min)

  surf0 <- pred_surv(fit = lasso_fit0,
                       S0 = S0,
                       x = data.test$X,
                       times = times,
                       lambda = lasso_fit0$lambda.min)

  pred_T_lasso <- surf1 - surf0
  pred_T_lasso
}
