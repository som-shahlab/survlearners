#' @title S-learner of lasso
#'
#' @description  S-learner, implemented via glmnet (lasso)
#'
#' @param X The baseline covariates
#' @param W The treatment variable (0 or 1)
#' @param Y The follow-up time
#' @param D The event indicator
#' @param times The prediction time of interest
#' @param alpha Mix tuning parameter for the elastic net
#' @param k.folds Number of folds for cross validation
#' @param foldid User-supplied foldid. Must have length equal to length(W). If provided, it overrides the k.folds option.
#' @param lambda User-supplied lambda sequence for cross validation
#' @param lambda.choice How to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty.factor User-supplied penalty factor, must be of length the same as number of features in X
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' times <- 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[ ,1] + (-0.5 - 1 * X[ ,2]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.sl.lasso.fit <- surv_sl_lasso(X, Y, W, D, times)
#' cate <- predict(surv.sl.lasso.fit)
#' cate.test <- predict(surv.sl.lasso.fit, X.test)
#' }
#' @return a surv_sl_lasso object
#' @export
surv_sl_lasso <- function(X, Y, W, D, times,
                          alpha = 1,
                          k.folds = NULL,
                          foldid = NULL,
                          lambda = NULL,
                          lambda.choice = "lambda.min",
                          penalty.factor = NULL) {

  input <- sanitize_input(X, Y, W, D)
  X <- input$X
  W <- input$W
  Y <- input$Y
  D <- input$D

  x.scl <- scale(X, center = TRUE, scale = TRUE)
  x.scl <- x.scl[ ,!is.na(colSums(x.scl)), drop = FALSE]

  nobs <- nrow(x.scl)
  pobs <- ncol(x.scl)

  if (is.null(foldid) || length(foldid) != length(W)) {

    if (!is.null(foldid) && length(foldid) != length(W)) {
      warning("supplied foldid does not have the same length ")
    }

    if (is.null(k.folds)) {
      k.folds <- floor(max(3, min(10, length(W) / 4)))
    }

    # fold ID for cross-validation; balance treatment assignments
    foldid <- sample(rep(seq(k.folds), length = length(W)))
  }

  x.scl.tilde <- cbind(as.numeric(2 * W - 1) * cbind(1, x.scl), x.scl)
  x.scl.pred1 <- cbind(1, x.scl, x.scl)
  x.scl.pred0 <- cbind(0, 0 * x.scl, x.scl)

  if (is.null(penalty.factor) || (length(penalty.factor) != pobs)) {
    if (!is.null(penalty.factor) && length(penalty.factor) != 2 * pobs + 1) {
      warning("penalty.factor supplied is not 1 plus 2 times the number of columns in X. Using all ones instead.")
    }
    penalty.factor <- c(0, rep(1, 2 * pobs))
  }
  x.scl.tilde <- as.matrix(data.frame(x.scl.tilde))
  tau.fit <- glmnet::cv.glmnet(x.scl.tilde,
                               survival::Surv(Y, D),
                               family = "cox",
                               foldid = foldid,
                               lambda = lambda,
                               penalty.factor = penalty.factor,
                               standardize = FALSE,
                               alpha = alpha)

  s.beta <- t(as.vector(coef(tau.fit, s = lambda.choice)))
  s.beta.adj <- c(0.5 * s.beta[1:(1 + dim(X)[2])], s.beta[(2 + dim(X)[2]):dim(x.scl.tilde)[2]])

  link1 <- exp(x.scl.pred1 %*% s.beta.adj)
  link0 <- exp(x.scl.pred0 %*% s.beta.adj)
  S0.t <- base_surv(fit = tau.fit,
                    Y = Y,
                    D = D,
                    X = x.scl.tilde,
                    lambda = tau.fit$lambda.min)
  index <- findInterval(times, S0.t$time)
  S0 <- S0.t[index, ]$survival
  surv1 <- S0 ^ exp(link1)
  surv0 <- S0 ^ exp(link0)

  tau.hat <- as.numeric(surv1 - surv0)

  ret <- list(tau.fit = tau.fit,
              x.org = x.scl.tilde,
              y.org = Y,
              D.org = D,
              beta.org = s.beta,
              s.beta = s.beta.adj,
              S0.t = S0.t,
              times = times,
              tau.hat = tau.hat,
              lambda.choice = lambda.choice)

  class(ret) <- "surv_sl_lasso"
  ret
}

#' predict for surv_sl_lasso
#'
#' get estimated tau(X) using the trained surv_sl_lasso model
#'
#' @param object A surv_sl_lasso object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param times The prediction time of interest
#' @param ... Additional arguments (currently not used)
#'
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' times <- 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[ ,1] + (-0.5 - 1 * X[ ,2]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.sl.lasso.fit <- surv_sl_lasso(X, Y, W, D, times)
#' cate <- predict(surv.sl.lasso.fit)
#' cate.test <- predict(surv.sl.lasso.fit, X.test)
#' }
#'
#' @return A vector of estimated conditional average treatment effects
#' @export
predict.surv_sl_lasso <- function(object,
                                  newdata = NULL,
                                  times = NULL,
                                  ...) {
  if (!is.null(newdata)) {
    newdata <- sanitize_x(newdata)
    newdata.scl <- scale(newdata, center = TRUE, scale = TRUE)
    newdata.scl <- newdata.scl[ ,!is.na(colSums(newdata.scl)), drop = FALSE]
    newdata.scl.pred1 <- cbind(1, newdata.scl, newdata.scl)
    newdata.scl.pred0 <- cbind(0, 0 * newdata.scl, newdata.scl)

    link1 <- exp(newdata.scl.pred1 %*% object$s.beta)
    link0 <- exp(newdata.scl.pred0 %*% object$s.beta)

    if (is.null(times)) {
    times <- object$times
    }
    index <- findInterval(times, object$S0.t$time)
    S0 <- object$S0.t[index, ]$survival

    surv1 <- S0^exp(link1)
    surv0 <- S0^exp(link0)

    tau.hat <- as.numeric(surv1 - surv0)
  }
  else {
    tau.hat <- object$tau.hat
  }
  return(tau.hat)
}
