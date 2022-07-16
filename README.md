[![ci](https://github.com/som-shahlab/survlearners/actions/workflows/main.yml/badge.svg)](https://github.com/som-shahlab/survlearners/actions/workflows/main.yml)

# Benchmarking Metalearners for Survival Data

This repository contains the extended implementations of five state-of-the-art metalearners for survival outcomes. Metalearners are specific meta-algortihms that leverage predictive models, e.g., machine learning models, to solve the causal task of estimating treatment heterogeneity. While metalearners (S-, T-, X-, M-, and R-learners) have been developed for uncensored continuous or binary data, suvrlearners provide their adapted executions to the survival setting through inverse probability of censoring weighting when combined with two popular machine learning strategies (Lasso and generalized random forests).

For a tutorial of these metalearners with mathematical underpinning, please refer to the Chapter named [Treatment heterogeneity for survival outcomes]() in the **Handbook of Matching and Weighting Adjustments for Causal Inference**. In addition, this chapter includes a benchmarking study of the five metaleaners via a comprehensive simulation study to assess the impacts of several important factors in the data generating process (the complexity of the baseline risk function, the complexity of the CATE function, the magnitude of the heterogeneity in treatment effects, the censoring mechanism, and the imbalance of treated and control units) on estimation performance. The summarized recommendations and considerations for choosing and applying metalearners and machine learning models for HTE estimation are incorporated into survlearners.

Some links for getting started

* The [tutorial](https://som-shahlab.github.io/survlearners/) contains usage examples and relevant references
* The replication code for the metalearner benchmarking study is available in [experiments](https://github.com/som-shahlab/survlearners/tree/master/experiments)

### Installation
The current version of survlearners can be installed from source using [remotes](https://cran.r-project.org/web/packages/remotes/index.html).

```R
remotes::install_github("som-shahlab/survlearners", subdir = "r-package/survlearners")
```
### Usage Examples

The following script demonstrates the application of R-learner with random survival forest and Lasso for estimating the heterogenous treatment effects on the absolute  difference scale, i.e., difference in survival probabilities at the median follow-up time. For examples of applying other types of metalearners, please refer to the R [documentation](https://som-shahlab.github.io/survlearners/reference/index.html).

```R
library(survlearners)

# Generate data
n <- 1000; p <- 25
t0 <- 0.2
Y.max <- 2
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 0.5)
numeratorT <- -log(runif(n))
T <- (numeratorT / exp(1 * X[ ,1, drop = FALSE] + (-0.5 - 1 * X[ ,2, drop = FALSE]) * W)) ^ 2
failure.time <- pmin(T, Y.max)
numeratorC <- -log(runif(n))
censor.time <- (numeratorC / (4 ^ 2)) ^ (1 / 2)
Y <- pmin(failure.time, censor.time)
D <- as.integer(failure.time <= censor.time)
n.test <- 500
X.test <- matrix(rnorm(n.test * p), n.test, p)

# Conditional average treatment effects (CATEs)
# Train a CATE model using R-learner with random survival forest (RSF) and Lasso
# RSF is used for estimating the nuisance parameters, e.g., mean outcome and censoring weights
# Lasso is used for estimating the target parameter CATE
fit <- surv_rl_grf_lasso(X, Y, W, D, t0, W.hat = 0.5, cen.fit = "survival.forest")

# Estimating conditional average treatment effects for the test sample
cate <- predict(fit, X.test)
```

### Funding

This work was supported by R01 HL144555 from the National Heart, Lung, and Blood Institute (NHLBI).

### References
Xu, Y., Ignatiadis, N., Sverdrup E., Fleming S., Wager, S.,  and Shah, N. (2022). **Handbook of Matching and Weighting Adjustments for Causal Inference. Chapter: Treatment Heterogeneity with Survival Outcomes.** Chapman \& Hall/CRC Press (forthcoming). [arXiv]()

Athey S., Tibshirani J., and Wager S. (2019). **Generalized random forests.** The Annals of Statistics, 47(2):1148–1178. [Paper](https://projecteuclid.org/journals/annals-of-statistics/volume-47/issue-2/Generalized-random-forests/10.1214/18-AOS1709.full)

Tibshirani R. (1997). **The lasso method for variable selection in the Cox model.** Statistics in medicine, 16(4):385–395. [Paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1097-0258(19970228)16:4%3C385::AID-SIM380%3E3.0.CO;2-3)

Van der Laan M. J. and Robins J. M. (2003). **Unified methods for censored longitudinal data and causality,** volume 5. Springer. [Paper](https://link.springer.com/book/10.1007/978-0-387-21700-0)

Tsiatis. A. A. (2006). Semiparametric theory and missing data. [Paper](https://link.springer.com/book/10.1007/0-387-37345-4)

Tian L., Alizadeh A., Gentles A., and Tibshirani R. (2014). **A simple method for estimating interactions between a treatment and a large number of covariates.** Journal of the American Statistical Association, 109. [Paper](https://www.tandfonline.com/doi/full/10.1080/01621459.2014.951443)

Kunzel S. R., Sekhon J. S., Bickel P. J., and Yu B. (2019). **Metalearners for estimating heterogeneous treatment effects using machine learning.** Proceedings of the national academy of sciences, 116(10):4156–4165, 2019. [paper](https://www.pnas.org/doi/abs/10.1073/pnas.1804597116)

Nie X. and Wager S. (2021). **Quasi-oracle estimation of heterogeneous treatment effects.** Biometrika, 108(2):299–319. [Paper](https://academic.oup.com/biomet/article/108/2/299/5911092?login=true)

Cui Y., Kosorok M. R., Sverdrup E., Wager S., and Zhu R. (2020). **Estimating heterogeneous treatment effects with right-censored data via causal survival forests.** [arXiv](https://arxiv.org/abs/2001.09887)
