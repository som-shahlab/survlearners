base_surv <- function(fit, Y, D, x, lambda){
  data <- data.frame(t_event=Y, event=D, x)
  tab <- data.frame(table(data[data$event == 1, "t_event"]))
  y <- as.numeric(as.character(sort(unique(tab[,1]))))
  d <- tab[,2]  # number of events at each unique time

  betaHat <- as.vector(fit$glmnet.fit$beta[,fit$lambda==lambda])
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

sanitize_x = function(x){
	# make sure x is a numeric matrix with named columns (for caret)
	if (!is.matrix(x) | !is.numeric(x) | any(is.na(x))) {
		stop("x must be a numeric matrix with no missing values")
	}
	colnames(x) = paste0("covariate_", 1:ncol(x))
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
