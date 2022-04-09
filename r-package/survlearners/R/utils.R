
#' @title Compute baseline survival probability from lasso models
#'
#' @description Compute baseline survival probability from lasso models under a Cox PH distribution
#'
#' @param fit A cv.glmnet object
#' @param Y The follow-up time
#' @param D The event indicator
#' @param x A marix of covariates
#' @param lambda The mixing tuning parameter in glmnet
#'
#' @return A two-column matrix vector of baseline survival probabilities and corresponding failure times.
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
#' data <- list(X = X, W = W, Y = Y, D = D)
#' data.test <- list(X = X, W = W, Y = Y, D = D)
#'
#' foldid <- sample(rep(seq(10), length = length(data$Y)))
#' lasso_fit <- glmnet::cv.glmnet(data$X,
#'                                 Surv(data$Y, data$D),
#'                                 family = "cox",
#'                                 alpha = 1,
#'                                 foldid = foldid)
#'
#' S0 <- base_surv(fit = lasso_fit,
#'                   Y = data$Y,
#'                   D = data$D,
#'                   x = data$X,
#'                   lambda = lasso_fit$lambda.min)
#' }
#' @export
base_surv <- function(fit, Y, D, x, lambda){
  data <- data.frame(t_event=Y, event=D, x)
  tab <- data.frame(table(data[data$event == 1, "t_event"]))
  y <- as.numeric(as.character(sort(unique(tab[,1]))))
  d <- tab[,2]  # number of events at each unique time

  betaHat <- as.vector((fit$glmnet.fit$beta)[,fit$lambda==lambda])
  h0 <- rep(NA, length(y))
  for(l in 1:length(y)){
    h0[l] <- d[l] / sum(exp(x[data$t_event >= y[l], rownames(fit$glmnet.fit$beta)] %*% betaHat))
  }

  S0 <- exp(-cumsum(h0))
  outcome <- data.frame(time=y,survival=S0)
  outcome
}

#' @title Compute survival probability from lasso models
#'
#' @description Compute survival probability from lasso models under a Cox PH distribution
#'
#' @param fit A cv.glmnet object
#' @param S0 The baseline hazard
#' @param x A marix of covariates
#' @param times The time of interest
#' @param lambda The mixing tuning parameter in glmnet
#'
#' @return A vector of survival probabilities at time t
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
#' data <- list(X = X, W = W, Y = Y, D = D)
#' data.test <- list(X = X, W = W, Y = Y, D = D)
#'
#' foldid <- sample(rep(seq(10), length = length(data$Y)))
#' lasso_fit <- glmnet::cv.glmnet(data$X,
#'                                 Surv(data$Y, data$D),
#'                                 family = "cox",
#'                                 alpha = 1,
#'                                 foldid = foldid)
#'
#' S0 <- base_surv(fit = lasso_fit,
#'                   Y = data$Y,
#'                   D = data$D,
#'                   x = data$X,
#'                   lambda = lasso_fit$lambda.min)
#'
#' surf <- pred_surv(fit = lasso_fit,
#'                     S0 = S0,
#'                      x = data.test$X,
#'                      times = times,
#'                      lambda = lasso_fit$lambda.min)
#' }
#' @export
pred_surv <- function(fit, S0, x, times, lambda){
  link <- predict(fit$glmnet.fit,x,type = "link")[,fit$lambda==lambda]
  colnames(link) <- NULL

  if(length(times)>1){
    S0_t <- rep(NA, length(times))
    for (i in 1:length(times)){
      S0_t[i] <- S0$survival[S0$time>=times[i]][1]
    }
  }else{
    S0_t <- S0$survival[S0$time>=times][1]
  }

  surv <- S0_t^exp(link)
  surv
}

pred_surv_preval <- function(fit, S0, times, lambda){
  link <- fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == lambda]
  colnames(link) <- NULL

  if(length(times)>1){
    S0_t <- rep(NA, length(times))
    for (i in 1:length(times)){
      S0_t[i] <- S0$survival[S0$time>=times[i]][1]
    }
  }else{
    S0_t <- S0$survival[S0$time>=times][1]
  }

  surv <- S0_t^exp(link)
  surv
}

#' @title Clean design matrix
#'
#' @description Make sure the covariate matrix is numeric and with no missing values
#'
#' @param x  A raw covariate matrix
#'
#' @return A cleaned analysis-ready covariate matrix
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' cleanX <- sanitize_x(X)
#' }
#' @export
sanitize_x = function(x){
	# make sure x is a numeric matrix with named columns (for caret)
	if (!is.matrix(x) | !is.numeric(x) | any(is.na(x))) {
		stop("x must be a numeric matrix with no missing values")
	}
	colnames(x) = paste0("covariate_", 1:ncol(x))
	return(x)
}

#' @title Clean input data
#'
#' @description Make sure the input covariate matrix and outcomes are numeric and with no missing values
#'
#' @param x  A raw covariate matrix
#' @param w  The treatment variable
#' @param y  The follow-up time
#' @param D  The event indicator
#'
#' @return A cleaned analysis-ready data set
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
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
#' data = sanitize_input(X, W, Y, D)
#' }
#' @export
sanitize_input = function(x,w,y,D) {
  x = sanitize_x(x)

  if (!is.numeric(w)) {
		stop("the input w should be a numeric vector")
  }
	if (is.numeric(w) & all(w %in% c(0,1))) {
		w = w==1
	}

	# make sure y is a numeric vector
	if (!is.numeric(y)) {
		stop("y should be a numeric vector")
	}

  if (!is.numeric(D)) {
    stop("the input D should be a numeric vector")
  }
  if (is.numeric(D) & all(D %in% c(0,1))) {
    D = D==1
  }

	# make sure the dimensions align
	if (length(y)!=nrow(x) | length(w)!=nrow(x) | length(D)!=nrow(x)) {
		stop("nrow(x), length(w), length(y), and length(D) should all be equal")
	}

	return(list(x=x,
	            w=w,
	            y=y,
	            D=D))
}
