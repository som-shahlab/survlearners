estimators <- list(
  surv_fl_grf = surv_fl_grf,
  surv_fl_lasso = surv_fl_lasso,
  surv_rl_grf_lasso = surv_rl_grf_lasso,
  surv_rl_f = surv_rl_f,
  surv_rl_lasso = surv_rl_lasso,
  surv_sl_coxph = surv_sl_coxph,
  surv_sl_grf = surv_sl_grf,
  surv_sl_lasso = surv_sl_lasso,
  surv_tl_coxph = surv_tl_coxph,
  surv_tl_grf = surv_tl_grf,
  surv_tl_lasso = surv_tl_lasso,
  surv_xl_grf_lasso = surv_xl_grf_lasso,
  surv_xl_grf = surv_xl_grf,
  surv_xl_lasso = surv_xl_lasso
)

test_that("estimators are ~ invariant to flipping treatment indicator", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  X.test <- matrix(rnorm(251 * p), 251, p)
  W <- rbinom(n, 1, 0.5)
  ft <- 1 + rexp(n, 0.5)
  ct <- 1 + rexp(n, 0.1)
  Y <- round(pmin(ft, ct), 2)
  D <- as.integer(ft <= ct)

  t0 <- median(Y)
  W.hat <- 0.5
  for (i in 1:length(estimators)) {
    estimator <- estimators[[i]]
    set.seed(i)
    if ("W.hat" %in% names(formals(estimator))) {
      fit <- estimator(X, Y, W, D, t0, W.hat = W.hat)
    } else {
      fit <- estimator(X, Y, W, D, t0)
    }

    set.seed(i)
    if ("W.hat" %in% names(formals(estimator))) {
      fit.f <- estimator(X, Y, 1 - W, D, t0, W.hat = W.hat)
    } else {
      fit.f <- estimator(X, Y, 1 - W, D, t0)
    }

    if (is.list(predict(fit.f)) == TRUE) {
      expect_equal(unlist(predict(fit)), -1 * unlist(predict(fit.f)), tolerance = 0.1, label = names(estimators)[i])
      expect_equal(unlist(predict(fit, X.test)), -1 * unlist(predict(fit.f, X.test)), tolerance = 0.1, label = names(estimators)[i])
    } else {
      expect_equal(predict(fit), -1 * predict(fit.f), tolerance = 0.1, label = names(estimators)[i])
      expect_equal(predict(fit, X.test), -1 * predict(fit.f, X.test), tolerance = 0.1, label = names(estimators)[i])
    }
  }
})

test_that("estimators are ~ invariant to shifting Y by constant", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  X.test <- matrix(rnorm(251 * p), 251, p)
  W <- rbinom(n, 1, 0.5)
  ft <- 1 + rexp(n, 0.5)
  ct <- 1 + rexp(n, 0.1)
  Y <- round(pmin(ft, ct), 2)
  D <- as.integer(ft <= ct)

  t0 <- median(Y)
  W.hat <- 0.5
  for (i in 1:length(estimators)) {
    estimator <- estimators[[i]]
    set.seed(i)
    if ("W.hat" %in% names(formals(estimator))) {
      fit <- estimator(X, Y, W, D, t0, W.hat = W.hat)
    } else {
      fit <- estimator(X, Y, W, D, t0)
    }

    set.seed(i)
    if ("W.hat" %in% names(formals(estimator))) {
      fit.f <- estimator(X, Y + 100, W, D, t0 + 100, W.hat = W.hat)
    } else {
      fit.f <- estimator(X, Y + 100, W, D, t0 + 100)
    }
    if (is.list(predict(fit.f)) == TRUE) {
      expect_equal(unlist(predict(fit)), unlist(predict(fit.f)), tolerance = 0.05, label = names(estimators)[i])
      expect_equal(unlist(predict(fit, X.test)), unlist(predict(fit.f, X.test)), tolerance = 0.05, label = names(estimators)[i])
    } else {
      expect_equal(predict(fit), predict(fit.f), tolerance = 0.05, label = names(estimators)[i])
      expect_equal(predict(fit, X.test), predict(fit.f, X.test), tolerance = 0.05, label = names(estimators)[i])
    }
  }
})
