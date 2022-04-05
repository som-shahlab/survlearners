# dgps.R - convenience script for generating simulation data for grf.

#' Generate causal forest data
#'
#' The following DGPs are available for benchmarking purposes:
#' \itemize{
#'  \item "simple": tau = max(X1, 0), e = 0.4 + 0.2 * 1{X1 > 0}.
#'  \item "aw1": equation (27) of https://arxiv.org/pdf/1510.04342.pdf
#'  \item "aw2": equation (28) of https://arxiv.org/pdf/1510.04342.pdf
#'  \item "aw3": confounding is from "aw1" and tau is from "aw2"
#'  \item "aw3reverse": Same as aw3, but HTEs anticorrelated with baseline
#'  \item "ai1": "Setup 1" from section 6 of https://arxiv.org/pdf/1504.01132.pdf
#'  \item "ai2": "Setup 2" from section 6 of https://arxiv.org/pdf/1504.01132.pdf
#'  \item "kunzel": "Simulation 1" from A.1 in https://arxiv.org/pdf/1706.03461.pdf
#'  \item "nw1": "Setup A" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
#'  \item "nw2": "Setup B" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
#'  \item "nw3": "Setup C" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
#'  \item "nw4": "Setup D" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
#'}
#'
#' Each DGP is parameterized by
#' X: observables,
#' m: conditional mean of Y,
#' tau: treatment effect,
#' e: propensity scores,
#' V: conditional variance of Y.
#'
#' The following rescaled data is returned
#' m = m / sd(m) * sigma.m,
#' tau = tau / sd(tau) * sigma.tau,
#' V = V / mean(V) * sigma.noise^2,
#' W = rbinom(e),
#' Y = m + (W - e) * tau + sqrt(V) + rnorm(n).
#'
#' @param n The number of observations.
#' @param p The number of covariates (note: the minimum varies by DGP).
#' @param sigma.m The standard deviation of the unconditional mean of Y. Default is 1.
#' @param sigma.tau The standard deviation of the treatment effect. Default  is 0.1.
#' @param sigma.noise The conditional variance of Y. Default is 1.
#' @param dgp The kind of dgp. Default is "simple".
#'
#' @return A list consisting of:
#'  X, Y, W, tau, m, e, dgp.
#'
#' @examples
#' \donttest{
#' # Generate simple benchmark data
#' data <- generate_causal_data(100, 5, dgp = "simple")
#' # Generate data from Wager and Athey (2018)
#' data <- generate_causal_data(100, 5, dgp = "aw1")
#' data2 <- generate_causal_data(100, 5, dgp = "aw2")
#' }
#' @export
generate_causal_data <- function(n, p, sigma.m = 1, sigma.tau = 0.1, sigma.noise = 1,
                                 dgp = c("simple", "aw1", "aw2", "aw3", "aw3reverse",
                                         "ai1", "ai2", "kunzel", "nw1", "nw2", "nw3", "nw4")) {
  # To add an additonal DGP, fill in the template below and add an entry to `dgp` and `.minp`.
  .minp <- c(simple=3, aw1=2, aw2=2, aw3=1, aw3reverse=1,
             ai1=2, ai2=6, kunzel=2, nw1=5, nw2=5, nw3=3, nw4=5)
  dgp <- match.arg(dgp)
  minp <- .minp[dgp]
  if (p < minp) {
    msg <- paste0("Selected dgp ", dgp, " requires a minimum of ", minp, " variables.")
    stop(msg)
  }

  if (dgp == "kunzel") {
    if (!("MASS" %in% utils::installed.packages())) {
      msg <- paste0("Selected dgp ", dgp, " requires the MASS library.")
      stop(msg)
    }
  }

  # Create data
  if (dgp == "simple") {
    X <- matrix(rnorm(n * p), n, p)
    tau <- pmax(X[, 1], 0)
    e <- 0.4 + 0.2 * (X[, 1] > 0)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- X[, 2] + pmin(X[, 3], 0) + e * tau
    V <- 1
  } else if (dgp == "aw1") {
    # equation (27) of https://arxiv.org/pdf/1510.04342.pdf
    X <- matrix(runif(n * p, min = 0, max = 1), n, p)
    tau <- rep(0, n)  # Treatment effect is zero
    e <- (1 / 4) * (1 + dbeta(X[, 1], 2, 4))  # Confounding
    W <- rbinom(n = n, size = 1, prob = e)
    m <- 2 * X[, 1] - 1 + e * tau
    V <- 1
  } else if (dgp == "aw2") {
    # equation (28) of https://arxiv.org/pdf/1510.04342.pdf
    X <- matrix(runif(n * p), n, p)
    zeta1 <- 1 + 1 / (1 + exp(-20 * (X[, 1] - (1 / 3))))
    zeta2 <- 1 + 1 / (1 + exp(-20 * (X[, 2] - (1 / 3))))
    tau <- zeta1 * zeta2
    e <- rep(0.5, n)  # Randomized trial (no confounding)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- e * tau
    V <- 1
  } else if (dgp == "aw3") {
    # section 6.2 in https://arxiv.org/pdf/1610.01271.pdf
    # (confounding from aw1, tau from aw2)
    X <- matrix(runif(n * p), n, p)
    zeta1 <- 1 + 1 / (1 + exp(-20 * (X[, 1] - (1 / 3))))
    zeta2 <- 1 + 1 / (1 + exp(-20 * (X[, 2] - (1 / 3))))
    tau <- zeta1 * zeta2
    e <- (1 / 4) * (1 + dbeta(X[, 1], 2, 4))  # Confounding
    W <- rbinom(n = n, size = 1, prob = e)
    m <- 2 * X[, 1] - 1 + e * tau
    V <- 1
  } else if (dgp == "aw3reverse") {
    # Same as aw3, but HTEs anticorrelated with baseline
    X <- matrix(runif(n * p), n, p)
    zeta1 <- 1 + 1 / (1 + exp(20 * (X[, 1] - (1 / 3))))
    zeta2 <- 1 + 1 / (1 + exp(20 * (X[, 2] - (1 / 3))))
    tau <- zeta1 * zeta2
    e <- (1 / 4) * (1 + dbeta(X[, 1], 2, 4))  # Confounding
    W <- rbinom(n = n, size = 1, prob = e)
    m <- 2 * X[, 1] - 1 + e * tau
    V <- 1
  } else if (dgp == "ai1") {
    X <- matrix(rnorm(n, p), n, p)
    nu_x <- 0.5 * X[, 1] + X[, 2]
    tau <- 0.25 * X[, 1]
    e <- rep(0.5, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- nu_x + e * tau
    V <- 0.1^2
  } else if (dgp == "ai2") {
    X <- matrix(rnorm(n, p), n, p)
    nu_x <- 0.5 * X[, 1] + 0.5 * X[, 2] + X[, 3] + X[, 4] + X[, 5] + X[, 6]
    tau <- 0.5 * ((X[, 1] > 0) * X[, 1] + (X[, 2] > 0) * X[, 2])
    e <- rep(0.5, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- nu_x + e * tau
    V <- 0.1^2
  } else if (dgp == "kunzel") {
    # "Simulation 1" from A.1 in https://arxiv.org/pdf/1706.03461.pdf
    # Extremely unbalanced treatment assignment, easy treatment effect.
    X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = stats::toeplitz(0.5^seq(0, p - 1)))
    tau <- 8 * (X[, 2] > 0.1)
    beta <- runif(p, -5, 5)
    mu_0 <- X %*% beta + 5 * (X[, 1] > 0.5) + rnorm(n = n)
    mu_1 <- mu_0 + tau + rnorm(n = n)
    e <- rep(0.01, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- c(W * mu_1 + (1 - W) * mu_0 - (W - e) * tau)
    V <- 1
  } else if (dgp == "nw1") {
    # "Setup A" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Difficult nuisance components, easy treatment effect function.
    X <- matrix(runif(n * p), n, p)
    tau <- (X[, 1] + X[, 2]) / 2
    eta <- 0.1
    e <- pmax(eta, pmin(sin(pi * X[, 1] * X[, 2]), 1 - eta))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- sin(pi * X[, 1] * X[, 2]) + 2 * (X[, 3] - 0.5)^2 + X[, 4] + 0.5 * X[, 5] + e * tau
    V <- 1
  } else if (dgp == "nw2") {
    # "Setup B" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Randomized trial
    X <- matrix(rnorm(n * p), n, p)
    tau <- X[,1] + log(1 + exp(X[, 2]))
    e <- rep(0.5, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- pmax(0, X[, 1] + X[, 2], X[, 3]) + pmax(0, X[, 4] + X[, 5]) + e * tau
    V <- 1
  } else if (dgp == "nw3") {
    # "Setup C" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Easy propensity score, strong confounding, difficult baseline,
    # constant treatment effect
    X <- matrix(rnorm(n * p), n, p)
    tau <- rep(1, n)
    e <- 1 / (1 + exp(X[, 2] + X[, 3]))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- 2 * log(1 + exp(X[, 1] + X[, 2] + X[, 3])) + e * tau
    V <- 1
  } else if (dgp == "nw4") {
    # "Setup D" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Unrelated treatment and control arms
    # (No upside to learning them jointly)
    X <- matrix(rnorm(n * p), n, p)
    tau <- pmax(X[, 1] + X[, 2] + X[, 3], 0) - pmax(X[, 4] + X[, 5], 0)
    e <- 1 / (1 + exp(-X[, 1]) + exp(-X[, 2]))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- (pmax(X[, 1] + X[, 2] + X[, 3], 0) + pmax(X[, 4] + X[, 5], 0)) / 2 + e * tau
    V <- 1
  }

  # Scale and return data (rescale if `m` and `tau` is not constant, the NA check is for when n=1)
  if (!is.na(sd(m)) & !(sd(m) == 0)) {
    m <- m / sd(m) * sigma.m
  }
  if (!is.na(sd(tau)) & !(sd(tau) == 0)) {
    tau <- tau / sd(tau) * sigma.tau
  }
  V <- V / mean(V) * sigma.noise^2
  Y <- m + (W - e) * tau + sqrt(V) * rnorm(n)
  out <- list(X = X, Y = Y, W = W, tau = tau, m = m, e = e, dgp = dgp)

  out
}

#' Simulate causal survival data
#'
#' The following DGPs are available for benchmarking purposes, T is the failure time
#' and C the censoring time:
#' \itemize{
#'   \item "simple1": T = X1*eps + W, C ~ U(0, 2) where eps ~ Exp(1) and Y.max = 1.
#'   \item  "type1": T is drawn from an accelerated failure time model and C from a Cox model (scenario 1 in https://arxiv.org/abs/2001.09887)
#'   \item  "type2": T is drawn from a proportional hazard model and C from a accelerated failure time (scenario 2 in https://arxiv.org/abs/2001.09887)
#'   \item  "type3": T and C are drawn from a Poisson distribution  (scenario 3 in https://arxiv.org/abs/2001.09887)
#'   \item  "type4": T and C are drawn from a Poisson distribution  (scenario 4 in https://arxiv.org/abs/2001.09887)
#'   \item  "type5": is similar to "type2" but with censoring generated from an accelerated failure time model.
#' }
#' @param n The number of samples.
#' @param p The number of covariates.
#' @param Y.max The maximum follow-up time (optional).
#' @param X The covariates (optional).
#' @param n.mc The number of monte carlo draws to estimate the treatment effect with. Default is 10000.
#' @param dgp The type of DGP.
#'
#' @return A list with entries:
#'  `X`: the covariates, `Y`: the event times, `W`: the treatment indicator, `D`: the censoring indicator,
#'  `cate`: the treatment effect estimated by monte carlo, `cate.sign`: the true sign of the cate for ITR comparison,
#'  `dgp`: the dgp name, `Y.max`: the maximum follow-up time.
#'
#' @examples
#' \donttest{
#' # Generate data
#' n <- 1000
#' p <- 5
#' data <- generate_causal_survival_data(n, p)
#' # Get true CATE on a test set
#' X.test <- matrix(seq(0, 1, length.out = 5), 5, p)
#' cate.test <- generate_causal_survival_data(n, p, X = X.test)$cate
#' }
#'
#' @export

# Simulation from Recursively Imputed Survival Trees (Zhu and Kosorok 2012)
library(mvtnorm)
generate_survival_tree_data <- function(n, p, Y.max = NULL, X = NULL, n.mc = 10000, times = NULL,
                                          dgp = c("scenario1", "scenario2", "scenario3", "scenario4", "scenario5")) {
  dgp <- match.arg(dgp)
  if (!is.null(X)) {
    p <- NCOL(X)
    n <- NROW(X)
  }

  if (dgp == "scenario1") {
    # p = 25; times = 0.8
    if (is.null(Y.max)) {
      Y.max <- 4
    }
    if (is.null(X)) {
      rho <- 0.9
      sigma <- matrix(NA, p, p)
      for (i in 1:p){
        for (j in 1:p){
          sigma[i, j] <- rho^(abs(i-j))
        }
      }
      X <- pnorm(rmvnorm(n, mean = rep(0, p), sigma = sigma))
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    mu <- 0.1*(rowSums(X[, 11:20])) + (0.6 + 0.2*X[,1] + 0.2*X[, 11])*abs(W-1)
    failure.time <- pmin(rexp(n, mu), Y.max)
    censor.time <- rexp(n, mean(mu)/2)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time) # censoring rate 35%
    catesp <- rep(NA, n)
    eps <- rnorm(n.mc)
    for (i in 1:n) {
      ft1 <- pmin(rexp(n.mc, 0.1*(sum(X[i, 11:20]))), Y.max)
      ft0 <- pmin(rexp(n.mc, 0.1*(sum(X[i, 11:20])) + (0.6 + 0.2*X[i,1] + 0.2*X[i, 11])), Y.max)
      catesp[i] <- mean((ft1 > times) - (ft0 > times))
    }
    catesp.sign <- sign(catesp)
  } else if (dgp == "scenario2") {
    # p=10; times = 0.4
    if (is.null(Y.max)) {
      Y.max <- 6
    }

    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    mu <- sin(X[, 1]*pi) + 2*abs(X[,2]-0.5) + X[, 3]^3 + (0.5 + sin(X[, 1]*pi) + 2*abs(X[,2]-0.5) + X[, 3]^3)*abs(W-1)
    failure.time <- pmin(rexp(n, mu), Y.max)
    censor.time <- runif(n, 0, Y.max)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time) # censoring rate = 12%
    catesp <- rep(NA, n)
    eps <- rnorm(n.mc)
    for (i in 1:n) {
      ft1 <- pmin(rexp(n.mc, sin(X[i, 1]*pi) + 2*abs(X[i,2]-0.5) + X[i, 3]^3), Y.max)
      ft0 <- pmin(rexp(n.mc, sin(X[i, 1]*pi) + 2*abs(X[i,2]-0.5) + X[i, 3]^3 + (0.5 + sin(X[i, 1]*pi) + 2*abs(X[i,2]-0.5) + X[i, 3]^3)), Y.max)
      catesp[i] <- mean((ft1 > times) - (ft0 > times))
    }
    catesp.sign <- sign(catesp)
  } else if (dgp == "scenario3") {
    # p=25; times = 2.3
    if (is.null(Y.max)) {
      Y.max <- 10
    }

    if (is.null(X)) {
      rho <- 0.75
      sigma <- matrix(NA, p, p)
      for (i in 1:p){
        for (j in 1:p){
          sigma[i, j] <- rho^(abs(i-j))
        }
      }
      X <- rmvnorm(n, mean = rep(0, p), sigma = sigma)
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    shape <- 0.5 + 0.3*abs(rowSums(X[, 11:15])) + abs(0.5 + X[, 1] + X[, 11])*W
    failure.time <- pmin(rgamma(n, shape = shape, scale = 2), Y.max)
    censor.time <- runif(n, 0, 1.5*Y.max)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)  # censoring rate = 23%
    catesp <- rep(NA, n)
    eps <- rnorm(n.mc)
    for (i in 1:n) {
      ft0 <- pmin(rgamma(n.mc, shape = 0.5 + 0.3*abs(sum(X[i, 11:15])), scale = 2), Y.max)
      ft1 <- pmin(rgamma(n.mc, shape = 0.5 + 0.3*abs(sum(X[i, 11:15])) + abs(0.5 + X[i, 1] + X[i, 11]), scale = 2), Y.max)
      catesp[i] <- mean((ft1 > times) - (ft0 > times))
    }
    catesp.sign <- sign(catesp)
  }else if (dgp == "scenario4") {
    # p=25; times = 1.3
      if (is.null(Y.max)) {
        Y.max <- 4
      }

      if (is.null(X)) {
        rho <- 0.75
        sigma <- matrix(NA, p, p)
        for (i in 1:p){
          for (j in 1:p){
            sigma[i, j] <- rho^(abs(i-j))
          }
        }
        X <- rmvnorm(n, mean = rep(0, p), sigma = sigma)
      }
      e <- (1 + dbeta(X[, 1], 2, 4)) / 4
      W <- rbinom(n, 1, e)
      mu <- 0.1*abs(rowSums(X[, 1:5])) + 0.1*abs(rowSums(X[, 21:25])) + abs(0.1 + 0.5*X[, 1] + 0.5*X[, 21])*W
      failure.time <- pmin(rlnorm(n, mu), Y.max)
      censor.time <- rlnorm(n, mu+0.5)
      Y <- pmin(failure.time, censor.time)
      D <- as.integer(failure.time <= censor.time) # censoring rate = 32%
      catesp <- rep(NA, n)
      eps <- rnorm(n.mc)
      for (i in 1:n) {
        ft0 <- pmin(rlnorm(n.mc, 0.1*abs(rowSums(X[, 1:5])) + 0.1*abs(rowSums(X[, 21:25]))), Y.max)
        ft1 <- pmin(rlnorm(n.mc, 0.1*abs(rowSums(X[, 1:5])) + 0.1*abs(rowSums(X[, 21:25]))) + abs(0.1 + 0.5*X[, 1] + 0.5*X[, 21]), Y.max)
        catesp[i] <- mean((ft1 > times) - (ft0 > times))
      }
      catesp.sign <- sign(catesp)
  } else if (dgp == "scenario5") {
    # p=10; times = 0.5
    if (is.null(Y.max)) {
      Y.max <- 2
    }

    if (is.null(X)) {
      rho <- 0.2
      sigma <- matrix(NA, p, p)
      for (i in 1:p){
        for (j in 1:p){
          sigma[i, j] <- rho^(abs(i-j))
        }
      }
      X <- pnorm(rmvnorm(n, mean = rep(0, p), sigma = sigma))
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    mu <- exp(rowSums(X[, 1:3]) + (0.5 + X[, 1] + X[, 2])*abs(W-1))/(1+exp(rowSums(X[, 1:3]) + (0.5 + X[, 1] + X[, 2])*abs(W-1)))
    failure.time <- pmin(rexp(n, mu), Y.max)
    censor.time <- rbinom(n, 1, mu/2)*failure.time/3
    censor.time[censor.time==0] <- failure.time[censor.time==0] + 1  # check with Nikos on this
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)  # censoring rate = 42%
    catesp <- rep(NA, n)
    eps <- rnorm(n.mc)
    for (i in 1:n) {
      ft1 <- pmin(rexp(n, exp(rowSums(X[, 1:3]))/(1+exp(rowSums(X[, 1:3])))), Y.max)
      ft0 <- pmin(rexp(n, exp(rowSums(X[, 1:3]) + (0.5 + X[, 1] + X[, 2]))/(1+exp(rowSums(X[, 1:3]) + (0.5 + X[, 1] + X[, 2])))), Y.max)
      catesp[i] <- mean((ft1 > times) - (ft0 > times)) # times = 0.5
    }
    catesp.sign <- sign(catesp)
  }
    list(X = X, Y = Y, W = W, D = D, catesp = catesp, catesp.sign = catesp.sign, dgp = dgp, Y.max = Y.max)
}

# R-learner simulation (survival version)
generate_R_learner_survival_data <- function(n, p, Y.max = NULL, X = NULL, n.mc = 10000, times = NULL,
                                             dgp = c("RLsurv1", "RLsurv2", "RLsurv3", "RLsurv4")) {
  dgp <- match.arg(dgp)
  if (!is.null(X)) {
    p <- NCOL(X)
    n <- NROW(X)
  }

  if (dgp == "RLsurv1") {
    # p = 25; times = 3.5
    if (is.null(Y.max)) {
      Y.max <- 15
    }
    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    e <- pmax(pmin(sin(pi * X[, 1] * X[, 2]), 1 - 0.1), 0.1)
    W <- rbinom(n, 1, e)
    b0 <- sin(pi * X[,1] * X[,2]) + 2 * (X[, 3] - 0.5)^2 + X[, 4] + 0.5 * X[, 5]
    tau <- 2*(X[, 1] + X[, 2])
    failure.time <- pmin(exp(b0 + (W - 0.5) * tau + rnorm(n)), Y.max)
    censor.time <- exp(1.5*b0 + (W - 0.5) * tau/2 + rnorm(n))
    Y <- pmin(failure.time, censor.time); median(Y)
    D <- as.integer(failure.time <= censor.time); table(D) # event rate ~ 25%
    cen <- ifelse(D == 0 & Y < times, 0, 1); table(cen) # censoring rate ~ 13%
    catesp <- rep(NA, n)
    for (i in 1:n) {
      ft1 <- pmin(exp(b0[i] + (1 - 0.5) * tau[i] + rnorm(n.mc)), Y.max)
      ft0 <- pmin(exp(b0[i] + (0 - 0.5) * tau[i] + rnorm(n.mc)), Y.max)
      catesp[i] <- mean((ft1 > times) - (ft0 > times))
    }
    catesp.sign <- sign(catesp)
  } else if (dgp == "RLsurv2") {
    # p=25; times = 1
    if (is.null(Y.max)) {
      Y.max <- 10
    }
    if (is.null(X)) {
      X <- matrix(rnorm(n * p), n, p)
    }
    W <- rbinom(n, 1, 0.5)
    b0 <- pmax((X[,1] + X[,2]), X[, 3], 0) + pmax(X[, 4] + X[, 5], 0)
    tau <- X[, 1] + log(1 + exp(X[, 2]))
    numerator <- -log(runif(n))
    failure.time <- pmin(numerator/(0.1 * exp(b0 + (0.5 - W) * tau)), Y.max)
    censor.time <- numerator/(0.2 * exp(X[,1] + X[,2] + X[,4] + X[,5]))
    Y <- pmin(failure.time, censor.time); median(Y)
    D <- as.integer(failure.time <= censor.time); table(D) # event rate ~ 28%
    cen <- ifelse(D == 0 & Y < times, 0, 1); table(cen)    # censoring rate ~ 17%
    catesp <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    for (i in 1:n) {
      ft1 <- pmin(numerator / (0.1 * exp(b0[i] + (0.5 - 1) * tau[i])), Y.max)
      ft0 <- pmin(numerator / (0.1 * exp(b0[i] + (0.5 - 0) * tau[i])), Y.max)
      catesp[i] <- mean((ft1 > times) - (ft0 > times))
    }
    catesp.sign <- sign(catesp)
  } else if (dgp == "RLsurv3") {
    # p=25; times = 1
    if (is.null(Y.max)) {
      Y.max <- 5
    }
    if (is.null(X)) {
      X <- matrix(rnorm(n * p), n, p)
    }
    e <- 1/(1 + exp(X[, 2] + X[, 3]))
    W <- rbinom(n, 1, e)
    b0 <- 2 * log(1 + exp(X[,1] + X[,2] + X[, 3]))
    tau <- 1
    numerator <- -log(runif(n))
    failure.time <- pmin(numerator^(1/2)/(0.1 * exp(b0 + (0.5 - W) * tau)), Y.max)  # rho = 2, increasing hazard function
    censor.time <- numerator/(0.2 * exp((X[,1] + X[,2])))
    Y <- pmin(failure.time, censor.time); median(Y)
    D <- as.integer(failure.time <= censor.time); table(D) # event rate ~ 30%
    cen <- ifelse(D == 0 & Y < times, 0, 1); table(cen)    # censoring rate ~ 11%
    catesp <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    for (i in 1:n) {
      ft1 <- pmin(numerator^(1/2) / (0.1 * exp(b0[i] + (0.5 - 1) * tau)), Y.max)
      ft0 <- pmin(numerator^(1/2) / (0.1 * exp(b0[i] + (0.5 - 0) * tau)), Y.max)
      catesp[i] <- mean((ft1 > times) - (ft0 > times))
    }
    catesp.sign <- sign(catesp)
  }else if (dgp == "RLsurv4") {
    # p = 50; times = 1
    if (is.null(Y.max)) {
      Y.max <- 6
    }

    if (is.null(X)) {
      X <- matrix(rnorm(n * p), n, p)
    }
    e <- 1/(1 + exp(-X[, 1]) + exp(-X[, 2]))
    W <- rbinom(n, 1, e)
    b0 <- (pmax(X[, 1] + X[, 2] + X[, 3], 0) + pmax(X[, 4] + X[, 5], 0))/2
    tau <- pmax(X[, 1] + X[, 2] + X[, 3], 0) - pmax(X[, 4] + X[, 5], 0)
    failure.time <- pmin(exp(b0 + (W - 0.5) * tau + rnorm(n)), Y.max)
    censor.time <- exp(1.5*b0 + (W - 0.5) * tau/2 + rnorm(n))
    Y <- pmin(failure.time, censor.time); median(Y)
    D <- as.integer(failure.time <= censor.time); table(D) # censoring rate = 30%
    catesp <- rep(NA, n)
    for (i in 1:n) {
      ft1 <- pmin(exp(b0[i] + (1 - 0.5) * tau[i] + rnorm(n.mc)), Y.max)
      ft0 <- pmin(exp(b0[i] + (0 - 0.5) * tau[i] + rnorm(n.mc)), Y.max)
      catesp[i] <- mean((ft1 > times) - (ft0 > times))
    }
    catesp.sign <- sign(catesp)
  }
  list(X = X, Y = Y, W = W, D = D, catesp = catesp, catesp.sign = catesp.sign, dgp = dgp, Y.max = Y.max)
}

# SPRINT-based simulation by Scotty
generate_sprint_survival_data <- function(n, p, Y.max = NULL, X = NULL, n.mc = 10000, times = NULL,
                                          dgp = c("ate0", "fixHR5", "hteHR5", "fixHR10", "hteHR10",
                                                  "fixHR20", "hteHR20", "fixHR50", "hteHR50")) {
  if (dgp == "ate0") {
    # sprint simulator by Scotty, ATE = 0
    results_dict <- simulate_parametric_SPRINT_data(n = n,
                                                    e = 0.5,
                                                    fixed_treatment_effect_strength = 0, # reflect HR = 1
                                                    variable_treatment_effect_type = 'ascvd_correlated',
                                                    variable_treatment_effect_strength = 0,
                                                    use_fit_with_treatment_interactions = FALSE,
                                                    end_time=365.25*3,
                                                    event_rate_multiplier = 1.0)
    X <- as.matrix(results_dict$X)
    Y <- results_dict$Y
    W <- results_dict$W
    D <- results_dict$D
    Y.max <- NULL
    catesp <- results_dict$tau_sc_diff
    catesp.sign <- 2 * as.numeric(catesp > 0.007) - 1  # ATE is positive, methods that can only estimate ATE show perfect classification if use 0 as the cutoff
    cate <- NULL
    cate.sign <- NULL
  }else if (dgp == "fixHR5") {
    # Heterogeneity on HR = 0, 5% event rate
    results_dict <- simulate_parametric_SPRINT_data(n = n,
                                                    e = 0.5,
                                                    fixed_treatment_effect_strength = -0.3, # reflect HR = 0.75 in SPRINT (NEJM 2015)
                                                    variable_treatment_effect_type = 'ascvd_correlated',
                                                    variable_treatment_effect_strength = 0, # CATE does not depend on X via ASCVD risk
                                                    use_fit_with_treatment_interactions = FALSE, # HR is constant across subjects
                                                    end_time=365.25*3,
                                                    event_rate_multiplier = 1.0)
    X <- as.matrix(results_dict$X)
    Y <- results_dict$Y
    W <- results_dict$W
    D <- results_dict$D
    Y.max <- NULL
    catesp <- results_dict$tau_sc_diff
    catesp.sign <- 2 * as.numeric(catesp > 0.007) - 1
    cate <- NULL
    cate.sign <- NULL
  }else if (dgp == "hteHR5"){
    # Heterogeneity on HR is NOT 0, 5% event rate
    results_dict <- simulate_parametric_SPRINT_data(n = n,
                                                    e = 0.5,
                                                    fixed_treatment_effect_strength = -0.3, # reflect HR = 0.75 in SPRINT (NEJM 2015)
                                                    variable_treatment_effect_type = 'ascvd_correlated',
                                                    variable_treatment_effect_strength = 0,
                                                    use_fit_with_treatment_interactions = TRUE, # HR varies across subjects, directly depends on X
                                                    end_time = 365.25*3,
                                                    event_rate_multiplier = 1.3) # adjust this to ensure 5% event rate
    X <- as.matrix(results_dict$X)
    Y <- results_dict$Y
    W <- results_dict$W
    D <- results_dict$D
    Y.max <- NULL
    catesp <- results_dict$tau_sc_diff
    catesp.sign <- 2 * as.numeric(catesp > 0.007) - 1
    cate <- NULL
    cate.sign <- NULL
  }else if (dgp == "fixHR10") {
    # Heterogeneity on HR = 0, ~10% event rate
    results_dict <- simulate_parametric_SPRINT_data(n = n,
                                                    e = 0.5,
                                                    fixed_treatment_effect_strength = -0.3, # reflect HR = 0.75 in SPRINT (NEJM 2015)
                                                    variable_treatment_effect_type = 'ascvd_correlated',
                                                    variable_treatment_effect_strength = 0,
                                                    use_fit_with_treatment_interactions = FALSE,
                                                    end_time=365.25*3,
                                                    event_rate_multiplier = 2)  # higher event rate
    X <- as.matrix(results_dict$X)
    Y <- results_dict$Y
    W <- results_dict$W
    D <- results_dict$D
    Y.max <- NULL
    catesp <- results_dict$tau_sc_diff
    catesp.sign <- 2*as.numeric(catesp>0.007)-1
    cate <- NULL
    cate.sign <- NULL
  }else if (dgp == "hteHR10") {
    # Heterogeneity on HR is NOT 0, 10% event rate
    results_dict <- simulate_parametric_SPRINT_data(n = n,
                                                    e = 0.5,
                                                    fixed_treatment_effect_strength = -0.3, # reflect HR = 0.75 in SPRINT (NEJM 2015)
                                                    variable_treatment_effect_type = 'ascvd_correlated',
                                                    variable_treatment_effect_strength = 0,
                                                    use_fit_with_treatment_interactions = TRUE,
                                                    end_time=365.25*3,
                                                    event_rate_multiplier = 2.5)
    X <- as.matrix(results_dict$X)
    Y <- results_dict$Y
    W <- results_dict$W
    D <- results_dict$D
    Y.max <- NULL
    catesp <- results_dict$tau_sc_diff
    catesp.sign <- 2*as.numeric(catesp>0.007)-1
    cate <- NULL
    cate.sign <- NULL
  }else if (dgp == "fixHR20") {
    # Heterogeneity on HR = 0, ~20% event rate
    results_dict <- simulate_parametric_SPRINT_data(n = n,
                                                    e = 0.5,
                                                    fixed_treatment_effect_strength = -0.3, # reflect HR = 0.75 in SPRINT (NEJM 2015)
                                                    variable_treatment_effect_type = 'ascvd_correlated',
                                                    variable_treatment_effect_strength = 0,
                                                    use_fit_with_treatment_interactions = FALSE,
                                                    end_time=365.25*3,
                                                    event_rate_multiplier = 3.5)  # higher event rate
    X <- as.matrix(results_dict$X)
    Y <- results_dict$Y
    W <- results_dict$W
    D <- results_dict$D
    Y.max <- NULL
    catesp <- results_dict$tau_sc_diff
    catesp.sign <- 2*as.numeric(catesp>0.007)-1
    cate <- NULL
    cate.sign <- NULL
  }else if (dgp == "hteHR20") {
    # Heterogeneity on HR is NOT 0, 20% event rate
    results_dict <- simulate_parametric_SPRINT_data(n = n,
                                                    e = 0.5,
                                                    fixed_treatment_effect_strength = -0.3, # reflect HR = 0.75 in SPRINT (NEJM 2015)
                                                    variable_treatment_effect_type = 'ascvd_correlated',
                                                    variable_treatment_effect_strength = 0,
                                                    use_fit_with_treatment_interactions = TRUE,
                                                    end_time=365.25*3,
                                                    event_rate_multiplier = 5.2)
    X <- as.matrix(results_dict$X)
    Y <- results_dict$Y
    W <- results_dict$W
    D <- results_dict$D
    Y.max <- NULL
    catesp <- results_dict$tau_sc_diff
    catesp.sign <- 2*as.numeric(catesp>0.007)-1
    cate <- NULL
    cate.sign <- NULL
  }else if (dgp == "fixHR50") {
    # Heterogeneity on HR = 0, ~50% event rate
    results_dict <- simulate_parametric_SPRINT_data(n = n,
                                                    e = 0.5,
                                                    fixed_treatment_effect_strength = -0.3, # reflect HR = 0.75 in SPRINT (NEJM 2015)
                                                    variable_treatment_effect_type = 'ascvd_correlated',
                                                    variable_treatment_effect_strength = 0,
                                                    use_fit_with_treatment_interactions = FALSE,
                                                    end_time=365.25*3,
                                                    event_rate_multiplier = 10)  # higher event rate
    X <- as.matrix(results_dict$X)
    Y <- results_dict$Y
    W <- results_dict$W
    D <- results_dict$D
    Y.max <- NULL
    catesp <- results_dict$tau_sc_diff
    catesp.sign <- 2*as.numeric(catesp>0.007)-1
    cate <- NULL
    cate.sign <- NULL
  }else if (dgp == "hteHR50") {
    # Heterogeneity on HR is NOT 0, 50% event rate
    results_dict <- simulate_parametric_SPRINT_data(n = n,
                                                    e = 0.5,
                                                    fixed_treatment_effect_strength = -0.3, # reflect HR = 0.75 in SPRINT (NEJM 2015)
                                                    variable_treatment_effect_type = 'ascvd_correlated',
                                                    variable_treatment_effect_strength = 0,
                                                    use_fit_with_treatment_interactions = TRUE,
                                                    end_time=365.25*3,
                                                    event_rate_multiplier = 17)
    X <- as.matrix(results_dict$X)
    Y <- results_dict$Y
    W <- results_dict$W
    D <- results_dict$D
    Y.max <- NULL
    catesp <- results_dict$tau_sc_diff
    catesp.sign <- 2*as.numeric(catesp>0.007)-1
    cate <- NULL
    cate.sign <- NULL
  }

  list(X = X, Y = Y, W = W, D = D, catesp = catesp, catesp.sign = catesp.sign, dgp = dgp, Y.max = Y.max)
}
# causal survival forest paper (Cui 2022)
generate_cui_data <- function(n, p, p_b, p_i, Y.max = NULL, X = NULL, n.mc = 10000, times = NULL,
                              dgp = c("simple1", "type1", "type2", "type3", "type4")) {
  #.minp <- c(simple1 = 1, type1 = 5, type2 = 5, type3 = 5, type4 = 5, type5 = 5, type6 = 6, type7 = 7)
  dgp <- match.arg(dgp)
  #minp <- .minp[dgp]
  if (!is.null(X)) {
    p <- NCOL(X)
    n <- NROW(X)
  }
  # if (p < minp) {
  #   stop(paste("Selected dgp", dgp, "requires a minimum of", minp, "variables."))
  # }

  if (dgp == "simple1") {
    if (is.null(Y.max)) {
      Y.max <- 1
    }
    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    W <- rbinom(n, 1, 0.5)
    failure.time <- pmin(rexp(n) * X[, 1] + W, Y.max)
    censor.time <- 2 * runif(n)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    temp <- rexp(n.mc)
    cate <- rep(NA, n)
    for (i in 1:n) {
      cate[i] <- mean(pmin(temp * X[i, 1] + 1, Y.max) - pmin(temp * X[i, 1], Y.max))
    }
    cate.sign = rep(1, n)
  } else if (dgp == "type1") {
    # Type 1 from https://arxiv.org/abs/2001.09887 (Cox PH censor time)

    if (is.null(Y.max)) {
      Y.max <- 1.5
    }
    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    I1 <- X[,1] < 0.5
    ft <- exp(-1.85 - 0.8 * I1 + 0.7 * sqrt(X[, 2]) + 0.2 * X[, 3] +
                (0.7 - 0.4 * I1 - 0.4 * sqrt(X[, 2])) * W + rnorm(n))
    failure.time <- pmin(ft, Y.max)
    numerator <- -log(runif(n))
    denominator <- exp(-1.75 - 0.5 * sqrt(X[, 2]) + 0.2 * X[, 3]  + (1.15 + 0.5 * I1 - 0.3 * sqrt(X[, 2])) * W)
    censor.time <- (numerator / denominator)^(1/2)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    #cate <- rep(NA, n)
    catesp <- rep(NA, n)
    eps <- rnorm(n.mc)
    for (i in 1:n) {
      ft0 <- exp(-1.85 - 0.8 * I1[i] + 0.7 * sqrt(X[i, 2]) + 0.2 * X[i, 3] + eps)
      ft1 <- exp(-1.85 - 0.8 * I1[i] + 0.7 * sqrt(X[i, 2]) + 0.2 * X[i, 3] +
                   0.7 - 0.4 * I1[i] - 0.4 * sqrt(X[i, 2]) + eps)
      catesp[i] <- mean((pmin(ft1, Y.max) > times) - (pmin(ft0, Y.max) > times))
    }
    catesp.sign <- sign(catesp)

  } else if (dgp == "type2") {
    # Type 2 from https://arxiv.org/abs/2001.09887 (Cox PH failure time)
    if (is.null(Y.max)) {
      Y.max <- 2
    }
    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    numerator <- -log(runif(n))
    cox.ft <- (numerator / exp(X[,1] + (-0.5 + X[,2]) * W))^2; summary(cox.ft)
    failure.time <- pmin(cox.ft, Y.max)
    censor.time <- 3 * runif(n)
    Y <- pmin(failure.time, censor.time); median(Y)
    D <- as.integer(failure.time <= censor.time); table(D) # event rate ~ 17%
    cen <- ifelse(D == 0 & Y < times, 0, 1); table(cen)    # censoring rate ~ 3%
    catesp <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    for (i in 1:n) {
      cox.ft0 <- (numerator / exp(X[i, 1] + (-0.5 + X[i, 2]) * 0))^2
      cox.ft1 <- (numerator / exp(X[i, 1] + (-0.5 + X[i, 2]) * 1))^2
      catesp[i] <- mean((pmin(cox.ft1, Y.max) > times) - (pmin(cox.ft0, Y.max) > times))
    }
    catesp.sign <- sign(catesp)
  } else if (dgp == "type3") {
    # Type 3 from https://arxiv.org/abs/2001.09887 (Poisson)
    if (is.null(Y.max)) {
      Y.max <- 15
    }
    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    lambda.failure <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * W
    failure.time <- pmin(rpois(n, lambda = lambda.failure), Y.max)
    lambda.censor <- 12 + log(1 + exp(X[, 3]))
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    catesp <- rep(NA, n)
    lambda.failure.0 <- X[, 2]^2 + X[, 3] + 6
    lambda.failure.1 <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3)
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
      catesp[i] <- mean((pmin(ft1, Y.max) > times) - (pmin(ft0, Y.max) > times))
    }
    cate.sign <- sign(sqrt(X[, 1]) - 0.3)
    catesp.sign <- sign(catesp)
  } else if (dgp == "type4") {
    # Type 4 from https://arxiv.org/abs/2001.09887 (Poisson)
    if (is.null(Y.max)) {
      Y.max <- 3
    }
    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    e <- 1 / ((1 + exp(-X[, 1])) * (1 + exp(-X[, 2])))
    W <- rbinom(n, 1, e)
    lambda.failure <- X[,2] + X[, 3] + pmax(0, X[, 1] - 0.3) * W
    failure.time <- pmin(rpois(n, lambda = lambda.failure), Y.max)
    lambda.censor <- 1 + log(1 + exp(X[, 3]))
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    catesp <- rep(NA, n)
    lambda.failure.0 <- X[,2] + X[, 3]
    lambda.failure.1 <- X[,2] + X[, 3] + pmax(0, X[, 1] - 0.3)
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
      catesp[i] <- mean((pmin(ft1, Y.max) > times) - (pmin(ft0, Y.max) > times))
    }
    cate.sign <- sign(pmax(0, X[, 1] - 0.3))
    catesp.sign <- sign(catesp)
    # For X1 < 0.3 the cate is zero so both (0, 1) are optimal, and we can ignore this subset.
    cate.sign[X[, 1] < 0.3] <- NA
  }
  list(X = X, Y = Y, W = W, D = D, catesp = catesp, catesp.sign = catesp.sign, dgp = dgp, Y.max = Y.max)
}

# Metalearners-Benchmark dgps
generate_tutorial_survival_data <- function(n, p, p_b = NULL, p_i = NULL, f_b = NULL, f_i = NULL,
                                            n.mc = 10000, times = NULL, Y.max = NULL, pi = 0.5,
                                            beta = 1, gamma = 1, rho = 2, cen_scale = 4,
                                            cenM = "indX", dgp = "fcomplex") {
  dgp <- match.arg(dgp)

  if (dgp == "fcomplex") {
    if (is.null(Y.max)) {
      Y.max <- 2
    }
    X <- matrix(rnorm(n * p), n, p)
    indcator <- function(x){
      as.numeric(x > 0.5)
    }
    NLX <- data.frame(apply(X, 2, FUN = indcator))           # nonlinear/binary version of X
    names(NLX) <- paste0("V", 1:25)
    NLXs <- model.matrix(~ V1 + V2:V3 + V4:V5 + V6:V7 + V8:V9 + V10:V11 + V12:V13 +
                           V14:V15 + V16:V17 + V18:V19 + V20:V21 + V22:V23 + V24:V25 - 1, NLX)
    W <- rbinom(n, 1, pi)
    numerator <- -log(runif(n))
    if(f_b == "L" & f_i == "L"){
      if(p_b == 1 & p_i == 1){
        cox.ft <- (numerator / exp(beta * X[,1] + (-0.5 - gamma * X[,2]) * W))^2
      }else if(p_b == p & p_i == 1){
        betah <- rep(beta/sqrt(p_b), p_b)
        cox.ft <- (numerator / exp(X %*% betah + (-0.5 - gamma * X[,2]) * W))^2
      }else if(p_b == p & p_i == p){
        betah <- rep(beta/sqrt(p_b), p_b); gammah <- rep(gamma/sqrt(p_i), p_i)
        cox.ft <- (numerator / exp(X %*% betah + (-0.5 - X %*% gammah) * W))^2
      }
    }else if(f_b == "NL" & f_i == "L"){
      if(p_b == 1 & p_i == 1){
        cox.ft <- (numerator / exp(beta * as.numeric(X[,1] > 0.5) + (-0.5 - gamma * X[,2]) * W))^2
      }else if(p_b == p & p_i == 1){
        betah <- c(0.99, rep(0.33, (p_b-1)/2))
        cox.ft <- (numerator / exp(NLXs %*% betah + (-0.5 - gamma * X[,2]) * W))^2
      }else if(p_b == p & p_i == p){
        betah <- c(0.99, rep(0.33, (p_b-1)/2)); gammah <- rep(gamma/sqrt(p_i), p_i)
        cox.ft <- (numerator / exp(NLXs %*% betah + (-0.5 - X %*% gammah) * W))^2
      }
    }else if(f_b == "NL" & f_i == "NL"){
      if(p_b == 1 & p_i == 1){
        cox.ft <- (numerator / exp(beta * as.numeric(X[,1] > 0.5) + (-0.5 - gamma * as.numeric(X[,2] > 0.5)) * W))^2
      }else if(p_b == p & p_i == 1){
        betah <- c(0.99, rep(0.33, (p_b-1)/2))
        cox.ft <- (numerator / exp(NLXs %*% betah + (-0.5 - gamma * as.numeric(X[,2] > 0.5)) * W))^2
      }else if(p_b == p & p_i == p){
        betah <- c(0.99, rep(0.33, (p_b-1)/2)); gammah <- c(0.99, rep(0.33, (p_i-1)/2))
        cox.ft <- (numerator / exp(NLXs %*% betah + (-0.5 - NLXs %*% gammah) * W))^2
      }
    }
    failure.time <- pmin(cox.ft, Y.max); summary(failure.time)

    # varying censoring rate by changing cen_scale and rho
    numeratorC <- -log(runif(n))
    if(cenM == "dX"){
      cen_scale <- exp(0.5 + 2 * X[,1] + (1 + 2 * X[,2]) * W)
    }
    censor.time <- (numeratorC/(cen_scale^rho))^(1/rho)
    Y <- pmin(failure.time, censor.time); median(Y)
    D <- as.integer(failure.time <= censor.time); table(D); summary(Y[D==1])
    cen <- ifelse(D == 0 & Y < times, 0, 1); table(cen)         # censoring rate = 0.3 at times = 0.2
    event <- ifelse(D == 1 & Y < times, 1, 0); table(event)     # observed event rate = 0.3 at times = 0.2

    # generate true CATEs
    mu0sp <- mu1sp <- catesp <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    if(f_b == "L" & f_i == "L"){
      if(p_b == 1 & p_i == 1){
        for (i in 1:n) {
          cox.ft0 <- (numerator / exp(beta * X[i, 1] + (-0.5 - gamma * X[i, 2]) * 0))^2
          cox.ft1 <- (numerator / exp(beta * X[i, 1] + (-0.5 - gamma * X[i, 2]) * 1))^2
          mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
          mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
          catesp[i] <- mu1sp[i] - mu0sp[i]
          }
      }else if(p_b == p & p_i == 1){
        for (i in 1:n) {
          cox.ft0 <- (numerator / as.vector(exp(t(X[i,]) %*% betah + (-0.5 - gamma * X[i, 2]) * 0)))^2
          cox.ft1 <- (numerator / as.vector(exp(t(X[i,]) %*% betah + (-0.5 - gamma * X[i, 2]) * 1)))^2
          mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
          mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
          catesp[i] <- mu1sp[i] - mu0sp[i]
        }
      }else if(p_b == p & p_i == p){
        for (i in 1:n) {
          cox.ft0 <- (numerator / as.vector(exp(t(X[i,]) %*% betah + (-0.5 - t(X[i,]) %*% gammah) * 0)))^2
          cox.ft1 <- (numerator / as.vector(exp(t(X[i,]) %*% betah + (-0.5 - t(X[i,]) %*% gammah) * 1)))^2
          mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
          mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
          catesp[i] <- mu1sp[i] - mu0sp[i]
        }
      }
    }else if(f_b == "NL" & f_i == "L"){
      if(p_b == 1 & p_i == 1){
        for (i in 1:n) {
          cox.ft0 <- (numerator / exp(beta * as.numeric(X[i,1] > 0.5) + (-0.5 - gamma * X[i, 2]) * 0))^2
          cox.ft1 <- (numerator / exp(beta * as.numeric(X[i,1] > 0.5) + (-0.5 - gamma * X[i, 2]) * 1))^2
          mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
          mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
          catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        }else if(p_b == p & p_i == 1){
          for (i in 1:n) {
            cox.ft0 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - gamma * X[i, 2]) * 0)))^2
            cox.ft1 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - gamma * X[i, 2]) * 1)))^2
            mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
            mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
            catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        }else if (p_b == p & p_i == p){
          for (i in 1:n) {
            cox.ft0 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - t(X[i,]) %*% gammah) * 0)))^2
            cox.ft1 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - t(X[i,]) %*% gammah) * 1)))^2
            mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
            mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
            catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        }
      }else if(f_b == "NL" & f_i == "NL"){
        if(p_b == 1 & p_i == 1){
          for (i in 1:n) {
            cox.ft0 <- (numerator / exp(beta * as.numeric(X[i,1] > 0.5) + (-0.5 - gamma * as.numeric(X[i,2] > 0.5)) * 0))^2
            cox.ft1 <- (numerator / exp(beta * as.numeric(X[i,1] > 0.5) + (-0.5 - gamma * as.numeric(X[i,2] > 0.5)) * 1))^2
            mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
            mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
            catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        }else if(p_b == p & p_i == 1){
          for (i in 1:n) {
            cox.ft0 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - gamma * as.numeric(X[i,2] > 0.5)) * 0)))^2
            cox.ft1 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - gamma * as.numeric(X[i,2] > 0.5)) * 1)))^2
            mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
            mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
            catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        }else if(p_b == p & p_i == p){
          for (i in 1:n) {
            cox.ft0 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - t(NLXs[i, ]) %*% gammah) * 0)))^2
            cox.ft1 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - t(NLXs[i, ]) %*% gammah) * 1)))^2
            mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
            mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
            catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        }
      }
    sd(catesp)/sd(mu0sp)        # heterogeneity of CATE relative to variation in baseline
    catesp.sign <- sign(catesp)
  }
  list(X = X, Y = Y, W = W, D = D, catesp = catesp, catesp.sign = catesp.sign, dgp = dgp, Y.max = Y.max)
}
