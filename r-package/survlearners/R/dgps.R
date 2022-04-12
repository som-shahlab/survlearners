# dgps.R - convenience script for generating simulation data for grf.


#' @title Metalearners-Benchmark dgps
#'
#' @description  Simulation data sets for comparing metalearners on estimating CATE in survival outcomes
#'
#' @param n The sample size
#' @param p The number of covariates
#' @param p.b The number of variables in the main effect function
#' @param p.i The number of variables in the interaction term
#' @param f.b The function form of the main effects (linear / nonlinear)
#' @param f.i The function form of the treatment-covariate interactions (linear / nonlinear)
#' @param pi The propensity score
#' @param beta The coefficients of variables in the main effect
#' @param gamma The coefficients of variables in the interaction
#' @param rho The shape parameter in Weibull distribution for censoring time
#' @param cen.scale The scale parameter in Weibull distribution for censoring time
#' @param cenM The complexity of censoring mechanism (dependent / independent to covariates)
#' @param n.mc The number of monte carlo draws to estimate the treatment effect with. Default is 10000.
#' @param times The time of interest
#' @param Y.max The maximum failure time
#' @param dgp The type of DGP
#' @examples
#' \donttest{
#' n <- 1000; p <- 25
#' data <- generate_tutorial_survival_data(n, p, p.b = 1, p.i = 1, f.b = "L", f.i = "L", times = 0.2)
#' cate.true <- generate_tutorial_survival_data(n, p, p.b = 1, p.i = 1, f.b = "L", f.i = "L", times = 0.2)$catesp
#' }
#' @return A simulation data set
#' @export
generate_tutorial_survival_data <- function(n, p, p.b = NULL, p.i = NULL, f.b = NULL, f.i = NULL,
                                            pi = 0.5, beta = 1, gamma = 1, rho = 2, cen.scale = 4, cenM = "indX",
                                            n.mc = 10000, times = NULL, Y.max = NULL, dgp = "fcomplex") {

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
    if (f.b == "L" & f.i == "L") {
      if (p.b == 1 & p.i == 1) {
        cox.ft <- (numerator / exp(beta * X[ ,1] + (-0.5 - gamma * X[ ,2]) * W))^2
      } else if (p.b == p & p.i == 1) {
        betah <- rep(beta / sqrt(p.b), p.b)
        cox.ft <- (numerator / exp(X %*% betah + (-0.5 - gamma * X[ ,2]) * W))^2
      } else if (p.b == p & p.i == p) {
        betah <- rep(beta / sqrt(p.b), p.b); gammah <- rep(gamma / sqrt(p.i), p.i)
        cox.ft <- (numerator / exp(X %*% betah + (-0.5 - X %*% gammah) * W))^2
      }
    } else if (f.b == "NL" & f.i == "L") {
      if (p.b == 1 & p.i == 1) {
        cox.ft <- (numerator / exp(beta * as.numeric(X[ ,1] > 0.5) + (-0.5 - gamma * X[ ,2]) * W))^2
      } else if (p.b == p & p.i == 1) {
        betah <- c(0.99, rep(0.33, (p.b - 1) / 2))
        cox.ft <- (numerator / exp(NLXs %*% betah + (-0.5 - gamma * X[ ,2]) * W))^2
      } else if (p.b == p & p.i == p) {
        betah <- c(0.99, rep(0.33, (p.b - 1) / 2)); gammah <- rep(gamma / sqrt(p.i), p.i)
        cox.ft <- (numerator / exp(NLXs %*% betah + (-0.5 - X %*% gammah) * W))^2
      }
    } else if (f.b == "NL" & f.i == "NL") {
      if (p.b == 1 & p.i == 1) {
        cox.ft <- (numerator / exp(beta * as.numeric(X[ ,1] > 0.5) + (-0.5 - gamma * as.numeric(X[ ,2] > 0.5)) * W))^2
      } else if (p.b == p & p.i == 1) {
        betah <- c(0.99, rep(0.33, (p.b - 1) / 2))
        cox.ft <- (numerator / exp(NLXs %*% betah + (-0.5 - gamma * as.numeric(X[ ,2] > 0.5)) * W))^2
      } else if (p.b == p & p.i == p) {
        betah <- c(0.99, rep(0.33, (p.b - 1) / 2)); gammah <- c(0.99, rep(0.33, (p.i - 1) / 2))
        cox.ft <- (numerator / exp(NLXs %*% betah + (-0.5 - NLXs %*% gammah) * W))^2
      }
    }
    failure.time <- pmin(cox.ft, Y.max); summary(failure.time)

    # varying censoring rate by changing cen.scale and rho
    numeratorC <- -log(runif(n))
    if (cenM == "dX") {
      cen.scale <- exp(0.5 + 2 * X[ ,1] + (1 + 2 * X[ ,2]) * W)
    }
    censor.time <- (numeratorC / (cen.scale ^ rho)) ^ (1 / rho)
    Y <- pmin(failure.time, censor.time); median(Y)
    D <- as.integer(failure.time <= censor.time); table(D); summary(Y[D == 1])
    cen <- ifelse(D == 0 & Y < times, 0, 1); table(cen)         # censoring rate = 0.3 at times = 0.2
    event <- ifelse(D == 1 & Y < times, 1, 0); table(event)     # observed event rate = 0.3 at times = 0.2

    # generate true CATEs
    mu0sp <- mu1sp <- catesp <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    if (f.b == "L" & f.i == "L") {
      if (p.b == 1 & p.i == 1) {
        for (i in 1:n) {
          cox.ft0 <- (numerator / exp(beta * X[i, 1] + (-0.5 - gamma * X[i, 2]) * 0))^2
          cox.ft1 <- (numerator / exp(beta * X[i, 1] + (-0.5 - gamma * X[i, 2]) * 1))^2
          mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
          mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
          catesp[i] <- mu1sp[i] - mu0sp[i]
          }
      } else if (p.b == p & p.i == 1) {
        for (i in 1:n) {
          cox.ft0 <- (numerator / as.vector(exp(t(X[i, ]) %*% betah + (-0.5 - gamma * X[i, 2]) * 0)))^2
          cox.ft1 <- (numerator / as.vector(exp(t(X[i, ]) %*% betah + (-0.5 - gamma * X[i, 2]) * 1)))^2
          mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
          mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
          catesp[i] <- mu1sp[i] - mu0sp[i]
        }
      } else if (p.b == p & p.i == p) {
        for (i in 1:n) {
          cox.ft0 <- (numerator / as.vector(exp(t(X[i, ]) %*% betah + (-0.5 - t(X[i, ]) %*% gammah) * 0)))^2
          cox.ft1 <- (numerator / as.vector(exp(t(X[i, ]) %*% betah + (-0.5 - t(X[i, ]) %*% gammah) * 1)))^2
          mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
          mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
          catesp[i] <- mu1sp[i] - mu0sp[i]
        }
      }
    } else if (f.b == "NL" & f.i == "L") {
      if (p.b == 1 & p.i == 1) {
        for (i in 1:n) {
          cox.ft0 <- (numerator / exp(beta * as.numeric(X[i,1] > 0.5) + (-0.5 - gamma * X[i, 2]) * 0))^2
          cox.ft1 <- (numerator / exp(beta * as.numeric(X[i,1] > 0.5) + (-0.5 - gamma * X[i, 2]) * 1))^2
          mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
          mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
          catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        } else if (p.b == p & p.i == 1) {
          for (i in 1:n) {
            cox.ft0 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - gamma * X[i, 2]) * 0)))^2
            cox.ft1 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - gamma * X[i, 2]) * 1)))^2
            mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
            mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
            catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        } else if (p.b == p & p.i == p) {
          for (i in 1:n) {
            cox.ft0 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - t(X[i, ]) %*% gammah) * 0)))^2
            cox.ft1 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - t(X[i, ]) %*% gammah) * 1)))^2
            mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
            mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
            catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        }
      } else if (f.b == "NL" & f.i == "NL") {
        if (p.b == 1 & p.i == 1) {
          for (i in 1:n) {
            cox.ft0 <- (numerator / exp(beta * as.numeric(X[i,1] > 0.5) + (-0.5 - gamma * as.numeric(X[i,2] > 0.5)) * 0))^2
            cox.ft1 <- (numerator / exp(beta * as.numeric(X[i,1] > 0.5) + (-0.5 - gamma * as.numeric(X[i,2] > 0.5)) * 1))^2
            mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
            mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
            catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        } else if (p.b == p & p.i == 1) {
          for (i in 1:n) {
            cox.ft0 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - gamma * as.numeric(X[i,2] > 0.5)) * 0)))^2
            cox.ft1 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - gamma * as.numeric(X[i,2] > 0.5)) * 1)))^2
            mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
            mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
            catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        } else if (p.b == p & p.i == p) {
          for (i in 1:n) {
            cox.ft0 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - t(NLXs[i, ]) %*% gammah) * 0)))^2
            cox.ft1 <- (numerator / as.vector(exp(t(NLXs[i, ]) %*% betah + (-0.5 - t(NLXs[i, ]) %*% gammah) * 1)))^2
            mu0sp[i] <- mean(pmin(cox.ft0, Y.max) > times)
            mu1sp[i] <-  mean(pmin(cox.ft1, Y.max) > times)
            catesp[i] <- mu1sp[i] - mu0sp[i]
          }
        }
      }
    sd(catesp) / sd(mu0sp)        # heterogeneity of CATE relative to variation in baseline
    catesp.sign <- sign(catesp)
  }
  list(X = X, Y = Y, W = W, D = D, catesp = catesp, catesp.sign = catesp.sign, dgp = dgp, Y.max = Y.max)
}
