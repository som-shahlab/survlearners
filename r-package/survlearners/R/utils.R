# For computing baseline hazard in coxph
base_surv <- function(fit, Y, D, X, lambda) {
  data <- data.frame(t.event = Y, event = D, X)
  tab <- data.frame(table(data[data$event == 1, "t.event"]))
  Y <- as.numeric(as.character(sort(unique(tab[ ,1]))))
  D <- tab[ ,2]  # number of events at each unique time

  betaHat <- as.vector((fit$glmnet.fit$beta)[ ,fit$lambda == lambda])
  h0 <- rep(NA, length(Y))
  for (l in 1:length(Y)) {
    h0[l] <- D[l] / sum(exp(X[data$t.event >= Y[l], rownames(fit$glmnet.fit$beta)] %*% betaHat))
  }

  S0 <- exp(-cumsum(h0))
  outcome <- data.frame(time = Y,survival = S0)
  outcome
}
pred_surv <- function(fit, S0, X, times, lambda) {
  link <- predict(fit$glmnet.fit,X,type = "link")[ ,fit$lambda == lambda]
  colnames(link) <- NULL

  if (length(times)>1) {
    S0.t <- rep(NA, length(times))
    for (i in 1:length(times)) {
      S0.t[i] <- S0$survival[S0$time >= times[i]][1]
    }
  } else {
    S0.t <- S0$survival[S0$time >= times][1]
  }

  surv <- S0.t^exp(link)
  surv
}
pred_surv_preval <- function(fit, S0, times, lambda) {
  link <- fit$fit.preval[ ,!is.na(colSums(fit$fit.preval))][ , fit$lambda[!is.na(colSums(fit$fit.preval))] == lambda]
  colnames(link) <- NULL

  if (length(times)>1) {
    S0.t <- rep(NA, length(times))
    for (i in 1:length(times)) {
      S0.t[i] <- S0$survival[S0$time >= times[i]][1]
    }
  } else {
    S0.t <- S0$survival[S0$time >= times][1]
  }

  surv <- S0.t^exp(link)
  surv
}


sanitize_x <- function(X) {
	# make sure X is a numeric matrix with named columns (for caret)
	if (!is.matrix(X) | !is.numeric(X) | any(is.na(X))) {
		stop("X must be a numeric matrix with no missing values")
	}
	colnames(X) <- paste0("covariate_", 1:ncol(X))
	return(X)
}

sanitize_input <- function(X, Y, W, D) {
  X <- sanitize_x(X)

  if (!is.numeric(W)) {
		stop("the input W should be a numeric vector")
  }
	if (is.numeric(W) & all(W %in% c(0,1))) {
		W <- W == 1
	}

	# make sure Y is a numeric vector
	if (!is.numeric(Y)) {
		stop("Y should be a numeric vector")
	}

  if (!is.numeric(D)) {
    stop("the input D should be a numeric vector")
  }
  if (is.numeric(D) & all(D %in% c(0,1))) {
    D <- D == 1
  }

	# make sure the dimensions align
	if (length(Y) != nrow(X) | length(W) != nrow(X) | length(D) != nrow(X)) {
		stop("nrow(X), length(W), length(Y), and length(D) should all be equal")
	}

	return(list(X = X,
	            Y = Y,
	            W = W,
	            D = D))
}
