#' @title R-learner of grf
#'
#' @description  R-learner, implemented via the grf package
#'
#' @param X The baseline covariates
#' @param Y The follow-up time
#' @param W The treatment variable (0 or 1)
#' @param D The event indicator
#' @param t0 The prediction time of interest
#' @param k.folds Number of folds for cross validation
#' @param W.hat Propensity score
#' @param Y.hat Conditional mean outcome E(Y|X)
#' @param C.hat Censoring weights
#' @param args.grf.nuisance Input arguments for a grf model that estimates nuisance parameters
#' @param args.grf.tau Input arguments for a grf model that estimates CATE
#' @param cen.fit The choice of model fitting for censoring
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' t0 <- 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[ ,1, drop = FALSE] + (-0.5 - 1 * X[ ,2, drop = FALSE]) * W)) ^ 2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.rl.grf.fit <- surv_rl_grf(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.rl.grf.fit)
#' cate.test <- predict(surv.rl.grf.fit, X.test)
#' }
#' @return A surv_rl_grf_fit object
#' @export
surv_rl_grf <- function(X, Y, W, D,
                        t0 = NULL,
                        k.folds = NULL,
                        W.hat = NULL,
                        Y.hat = NULL,
                        C.hat = NULL,
                        args.grf.nuisance = list(),
                        args.grf.tau = list(),
                        cen.fit = "Kaplan-Meier") {

  input <- sanitize_input(X, Y, W, D)
  X <- input$X
  W <- as.numeric(input$W)
  Y <- input$Y
  D <- input$D
  nobs <- nrow(X)
  pobs <- ncol(X)

  if (is.null(k.folds)) {
    k.folds <- floor(max(3, min(10, length(Y) / 4)))
  }

  if (is.null(W.hat)) {
    stop("propensity score needs to be supplied")
  } else if (length(W.hat) == 1) {
    W.hat <- rep(W.hat, nrow(X))
  } else if (length(W.hat) != nrow(X)) {
    stop("W.hat has incorrect length.")
  }

  args.grf.nuisance <- list(failure.times = NULL,
                            num.trees = max(50, 2000 / 4),
                            sample.weights = NULL,
                            clusters = NULL,
                            equalize.cluster.weights = FALSE,
                            sample.fraction = 0.5,
                            mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                            min.node.size = 15,
                            honesty = TRUE,
                            honesty.fraction = 0.5,
                            honesty.prune.leaves = TRUE,
                            alpha = 0.05,
                            prediction.type = "Nelson-Aalen",
                            compute.oob.predictions = FALSE,
                            num.threads = NULL,
                            seed = runif(1, 0, .Machine$integer.max))

  if (is.null(Y.hat)) {
    y.fit <- do.call(grf::survival_forest, c(list(X = cbind(X, W), Y = Y, D = D), args.grf.nuisance))
    y.fit[["X.orig"]][ , ncol(X) + 1] <- rep(1, nrow(X))
    S1.hat <- predict(y.fit)$predictions
    y.fit[["X.orig"]][ , ncol(X) + 1] <- rep(0, nrow(X))
    S0.hat <- predict(y.fit)$predictions
    y.fit[["X.orig"]][ , ncol(X) + 1] <- W

    t0.index <- findInterval(t0, y.fit$failure.times)
    if (t0.index == 0) {
      Y.hat <- rep(1, nrow(X))
    } else {
      surf1 <- S1.hat[ , t0.index, drop = FALSE]
      surf0 <- S0.hat[ , t0.index, drop = FALSE]
      Y.hat  <- W.hat * surf1 + (1 - W.hat) * surf0
    }
  } else {
    y.fit <- NULL
  }

  if (is.null(C.hat)) {
    Q <- as.numeric(D == 1 | Y > t0)         # indicator for uncensored at t0
    U <- pmin(Y, t0)                         # truncated follow-up time by t0
    if (cen.fit == "Kaplan-Meier") {
      fold.id <- sample(rep(seq(k.folds), length = nrow(X)))
      C.hat <- rep(NA, length(fold.id))
      for (z in 1:k.folds) {
        c.fit <- survival::survfit(survival::Surv(Y[!fold.id == z], 1 - Q[!fold.id == z]) ~ 1)
        C.hat[fold.id == z] <- summary(c.fit, times = U[fold.id == z])$surv
      }
    } else if (cen.fit == "survival.forest") {
      args.grf.nuisance$compute.oob.predictions <- TRUE
      c.fit <- do.call(grf::survival_forest, c(list(X = cbind(W, X), Y = Y, D = 1 - Q), args.grf.nuisance))
      C.hat <- predict(c.fit)$predictions
      cen.times.index <- findInterval(U, c.fit$failure.times)
      C.hat <- C.hat[cbind(1:length(U), cen.times.index)]
    }
  } else {
    c.fit <- NULL
  }

  # CATE function
  D.t0 <- D
  D.t0[D == 1 & Y > t0] <- 0
  D.t0 <- D.t0[D == 1 | Y > t0]
  Y.hat.t0 <- Y.hat[D == 1 | Y > t0]
  W.t0 <- W[D == 1 | Y > t0]
  W.hat.t0 <- W.hat[D == 1 | Y > t0]
  X.t0 <- X[D == 1 | Y > t0,, drop = FALSE]
  C.hat.t0 <- C.hat[D == 1 | Y > t0]

  y.tilde <- (1 - D.t0) - Y.hat.t0
  w.tilde <- W.t0 - W.hat.t0
  pseudo.outcome <- y.tilde / w.tilde
  sample.weights <- w.tilde^2 / C.hat.t0

  args.grf.tau <- list(sample.weights = sample.weights,
                       num.trees = 2000,
                       clusters = NULL,
                       sample.fraction = 0.5,
                       mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                       min.node.size = 5,
                       honesty = TRUE,
                       honesty.fraction = 0.5,
                       honesty.prune.leaves = TRUE,
                       imbalance.penalty = 0,
                       ci.group.size = 2,
                       compute.oob.predictions = TRUE,
                       num.threads = NULL,
                       seed = runif(1, 0, .Machine$integer.max))

  tau.fit <- do.call(grf::regression_forest, c(list(X = X.t0, Y = pseudo.outcome), args.grf.tau))
  tau.hat <- predict(tau.fit, data.frame(X))

  ret <- list(tau.fit = tau.fit,
              pseudo.outcome = pseudo.outcome,
              sample.weights = sample.weights,
              y.fit = y.fit,
              c.fit = c.fit,
              W.hat = W.hat,
              Y.hat = Y.hat,
              C.hat = C.hat,
              tau.hat = tau.hat)
  class(ret) <- "surv_rl_grf"
  ret
}

#' predict for surv_rl_grf
#'
#' get estimated tau(X) using the trained surv_rl_grf model
#'
#' @param object A surv_rl_grf object
#' @param newdata Covariate matrix to make predictions on. If null, return the tau(X) predictions on the training data
#' @param ... Additional arguments (currently not used)
#'
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' t0 <- 0.2
#' Y.max <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' numeratorT <- -log(runif(n))
#' T <- (numeratorT / exp(1 * X[ ,1, drop = FALSE] + (-0.5 - 1 * X[ ,2, drop = FALSE]) * W))^2
#' failure.time <- pmin(T, Y.max)
#' numeratorC <- -log(runif(n))
#' censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' n.test <- 500
#' X.test <- matrix(rnorm(n.test * p), n.test, p)
#'
#' surv.rl.grf.fit <- surv_rl_grf(X, Y, W, D, t0, W.hat = 0.5)
#' cate <- predict(surv.rl.grf.fit)
#' cate.test <- predict(surv.rl.grf.fit, X.test)
#' }
#'
#' @return A vector of predicted conditional average treatment effects
#' @export
predict.surv_rl_grf <- function(object,
                                newdata = NULL,
                                ...) {
  if (!is.null(newdata)) {
    newdata <- sanitize_x(newdata)
    tau.hat <- predict(object$tau.fit, newdata)$predictions
  } else {
    tau.hat <- object$tau.hat
  }
  return(tau.hat)
}
