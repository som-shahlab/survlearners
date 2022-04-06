# For computing baseline hazard in coxph
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


# For thresholding propensity scores
trim = function(x, min, max) {
	x[x>max] = max
	x[x<min] = min
	return(x)
}

sanitize_x = function(x){
	# make sure x is a numeric matrix with named columns (for caret)
	if (!is.matrix(x) | !is.numeric(x) | any(is.na(x))) {
		stop("x must be a numeric matrix with no missing values")
	}
	colnames(x) = stringr::str_c("covariate_", 1:ncol(x))
	return(x)
}

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

#' @title data simulation
#'
#' @description Generates a dataset of size \eqn{n} that can be used to experiment with the learners and meta-learners
#' in this package.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' data = data_simulation(500) # draw a sample
#' @export
data_simulation = function(n) {
  x = stats::model.matrix(~.-1, data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n), "covariate_3" = rnorm(n), "covariate_4" = rnorm(n), "covariate_5" = rnorm(n), "covariate_6" = rnorm(n)))
	p = 0.5
	w = as.numeric(rbinom(n,1,p)==1)
  m = pmax(0, x[,1] + x[,2], x[,3]) + pmax(0, x[,4] + x[,5])
  tau = x[,1] + log(1 + exp(x[,2]))
	mu1 = m + tau/2
	mu0 = m - tau/2
	y = w*mu1 + (1-w) * mu0 + 0.5*rnorm(n)
	list(x=x, w=w, y=y, p=p, m=m, mu0=mu0, mu1=mu1, tau=tau)
}

#' @title Toy data simulation
#'
#' @description Generates a toy dataset of size \eqn{n} that can be used to experiment with the learners and meta-learners
#' in this package. The generative process should be easy to learn with linear methods.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' data = easy_toy_data_simulation(500) # draw a sample
#' @export
easy_toy_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n)))
	p = rep(0.5, n)
	w = as.numeric(rbinom(n,1,p)==1)
	tau = x %*% c(1,1)
	m = x %*% c(0.5,-0.5)
	mu1 = m + tau/2
	mu0 = m - tau/2
	y = (m + tau/2*(2*w-1))[,1]
	list(x=x, w=w, y=y, p=p, m=m, mu0=mu0, mu1=mu1, tau=tau)
}
#' @title Toy data simulation with continuous treatments
#'
#' @description Generates a toy dataset of size \eqn{n} that can be used to experiment with the learners and meta-learners
#' in this package. The generative process should be easy to learn with linear methods.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' data = continuous_toy_data_simulation(500) # draw a sample
#' @export
continuous_toy_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n)))
	w = runif(n, 0, 1)
	tau = x %*% c(1,1)
	m = x %*% c(0.5,-0.5)
	mu1 = m + tau/2
	mu0 = m - tau/2
	y = (m + tau/2*(2*w-1))[,1]
	list(x=x, w=w, y=y,  m=m, mu0=mu0, mu1=mu1, tau=tau)
}

#' @title Toy data simulation for T-learner
#'
#' @description Generates a toy dataset of size \eqn{n} that can be used to experiment with T-learners
#' in this package. The generative process should be easy to learn with linear methods.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' data = t_toy_data_simulation(500) # draw a sample
#' @export
t_toy_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n)))
	p = rep(0.5, n)
	w = as.numeric(rbinom(n,1,p)==1)
	mu1 = sin(x[,1] * 2)
	mu0 = x[,2] * 3 + 10
	y = w * mu1 + (1-w) * mu0
	tau = mu1 - mu0
	m = p * mu1 + (1-p) * mu0
	list(x=x, w=w, y=y, p=p, m=m, mu0=mu0, mu1=mu1, tau=tau)
}
#' @title data simulation for T-learner
#'
#' @description Generates a dataset of size \eqn{n} that can be used to experiment with T-learners
#' in this package. The generative process should be easy to learn with linear methods.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' data = t_data_simulation(500) # draw a sample
#' @export
t_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n), "covariate_3" = rnorm(n), "covariate_4" = rnorm(n), "covariate_5" = rnorm(n), "covariate_6" = rnorm(n)))
  p = 1/(1 + exp(-x[,1]) + exp(-x[,2]))
	w = as.numeric(rbinom(n,1,p)==1)
  b = (pmax(x[,1] + x[,2] + x[,3], 0) + pmax(x[,4] + x[,5], 0)) / 2
  tau = pmax(x[,1] + x[,2] + x[,3], 0) - pmax(x[,4] + x[,5], 0)

  mu1 = b + 0.5 * tau
  mu0 = b - 0.5 * tau

	y = w * mu1 + (1-w) * mu0 + rnorm(n)
	m = p * mu1 + (1-p) * mu0
	list(x=x, w=w, y=y, p=p, m=m, mu0=mu0, mu1=mu1, tau=tau)
}
#' @title helper function for result comparison
#' @description helper function to check if the learned tau_hat is within some bounded error from the groundtruth
#' @param tau_hat user-supplied treatment effect estimate
#' @param sim_data simulated groundtruth data
#' @param mse user-supplied error tolerance
#'
#' @export
meta_learner_tests = function(tau_hat, sim_data, mse=0.01) {
  expect_equal(length(tau_hat), length(sim_data$tau))
  expect_equal(is.numeric(tau_hat), TRUE)
  learner_mse = mean((tau_hat - sim_data$tau)^2)
  print(learner_mse)
  expect_equal(learner_mse<mse, TRUE)
}

#' @title helper function for testing the code runs
#' @description helper function to check if the learned tau_hat is numeric
#' @param tau_hat user-supplied treatment effect estimate
#' @param sim_data simulated groundtruth data
#' @param mse user-supplied error tolerance
#'
#' @export
simple_meta_learner_tests = function(tau_hat, sim_data, mse=0.01) {
  expect_equal(length(tau_hat), length(sim_data$tau))
  expect_equal(is.numeric(tau_hat), TRUE)
}

#' @title helper function for testing treatment effect is invariant when outcome adds 1
#' @description helper function to test treatment effect is invariant when outcome adds 1
#' @param tau_hat user-supplied treatment effect estimate
#' @param tau_hat_1 user-supplied treatment effect estimate for the setting with outcome added 1
#' @param mean_err error tolerance on the mean difference bewteen tau_hat and tau_hat_1
#'
#' @export
invariate_add_tests = function(tau_hat, tau_hat_1, mean_err=0.15) {
  print(abs(mean(tau_hat - tau_hat_1)))
  expect_equal(abs(mean(tau_hat - tau_hat_1)) < mean_err, TRUE)
}

#' @title helper function for testing treatment effect is invariant with a factor of 2 when outcome is multiplied with 2
#' @description helper function to test treatment effect is invariant with a factor of 2 when outcome is multiplied with 2
#' @param tau_hat user-supplied treatment effect estimate
#' @param tau_hat_2 user-supplied treatment effect estimate for the setting with outcome is multiplied with 2
#' @param mean_err error tolerance on the mean between 2x tau_hat and tau_hat_2
#'
#' @export
invariate_mult_tests = function(tau_hat, tau_hat_2, mean_err = 0.1) {
  print(abs(mean(2*tau_hat - tau_hat_2)) )
  expect_equal(abs(mean(2*tau_hat - tau_hat_2)) < mean_err, TRUE)
}
