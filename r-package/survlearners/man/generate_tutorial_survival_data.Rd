% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgps.R
\name{generate_tutorial_survival_data}
\alias{generate_tutorial_survival_data}
\title{Metalearners-Benchmark dgps}
\usage{
generate_tutorial_survival_data(
  n,
  p,
  p_b = NULL,
  p_i = NULL,
  f_b = NULL,
  f_i = NULL,
  pi = 0.5,
  beta = 1,
  gamma = 1,
  rho = 2,
  cen_scale = 4,
  cenM = "indX",
  n.mc = 10000,
  times = NULL,
  Y.max = NULL,
  dgp = "fcomplex"
)
}
\arguments{
\item{n}{The sample size}

\item{p}{The number of covariates}

\item{p_b}{The number of variables in the main effect function}

\item{p_i}{The propensity score}

\item{f_b}{The function form of the main effects (linear/nonlinear)}

\item{f_i}{The function form of the treatment-covariate interactions (linear/nonlinear)}

\item{beta}{The coefficients of variables in the main effect}

\item{gamma}{The coefficients of variables in the interaction}

\item{rho}{The shape parameter in Weibull distribution for censoring time}

\item{cen_scale}{The scale parameter in Weibull distribution for censoring time}

\item{cenM}{The complexity of censoring mechanism (dependent/independent to covariates)}

\item{n.mc}{The number of monte carlo draws to estimate the treatment effect with. Default is 10000.}

\item{times}{The time of interest}

\item{Y.max}{The maximum failure time}

\item{dgp}{The type of DGP}
}
\value{
A simulation data set
}
\description{
Simulation data sets for comparing metalearners on estimating CATE in survival outcomes
}
\examples{
\donttest{
n = 1000; p = 25
data = generate_tutorial_survival_data(n, p, p_b = 1, p_i = 1, f_b = "L", f_i = "L", times = 0.2)
cate.true = generate_tutorial_survival_data(n, p, p_b = 1, p_i = 1, f_b = "L", f_i = "L", times = 0.2)$catesp
}
}