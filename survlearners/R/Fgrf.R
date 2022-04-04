#'  Fit a pollinated transformed outcome (PTO) forest model
#'
#' @param x matrix of covariates
#' @param tx vector of treatment indicators (0 or 1)
#' @param y vector of response values
#' @param pscore vector of propensity scores
#' @param num.trees number of trees for transformed outcome forest
#' @param mtry number of variables to possibly split at in each node
#' @param min.node.size minimum node size for transformed outcome forest
#' @param postprocess logical: should optional post-processing random forest be
#'  fit at end?
#' @param verbose logical: should progress be printed to console?
#'
#' @return an object of class \code{PTOforest} with attributes:
#'  \itemize{
#'    \item x: matrix of covariates supplied by function call
#'    \item pscore: vector of propensity score supplied by function call
#'    \item postprocess: logical supplied by function call
#'    \item TOfit: fitted random forest on transformed outcomes
#'    \item PTOfit1: TOfit pollinated with treatment-arm outcomes
#'    \item PTOfit0: TOfit pollinated with control-arm outcomes
#'    \item postfit: post-processing random forest summarizing results
#'  }
#'
#' @examples
#'# Randomized experiment example
#'
#'n = 100 # number of training-set patients to simulate
#'p = 10  # number of features for each training-set patient
#'
#'# Simulate data
#'x = matrix(rnorm(n * p), nrow = n, ncol = p) # simulate covariate matrix
#'tx_effect = x[, 1] + (x[, 2] > 0) # simple heterogeneous treatment effect
#'tx = rbinom(n, size = 1, p = 0.5) # random treatment assignment
#'y = rowMeans(x) + tx * tx_effect + rnorm(n, sd = 0.001) # simulate response
#'
#'# Estimate PTO forest model
#'#fit_pto = PTOforest(x, tx, y)
#'#pred_pto = predict(fit_pto, newx = x)
#'
#'# Visualize results
#'#plot(tx_effect, pred_pto, main = 'PTO forest',
#'#  xlab = 'True treatment effect', ylab = 'Estimated treatment effect')
#'#abline(0, 1, lty = 2)
#'
#' @export


Fgrf = function(x, tx, y, pscore = rep(.5, nrow(x)),
                weight, num.trees = 2000, alpha = 0.05,
                meta_learner = TRUE, verbose = FALSE) {

  # Input sanitization

  x = as.matrix(x)

  if (nrow(x) != length(tx)) {
    stop('nrow(x) does not match length(tx)')

  } else if (nrow(x) != length(y)) {
    stop('nrow(x) does not match length(y)')

  } else if (!is.numeric(x)) {
    stop('x must be numeric matrix')

  } else if (!is.numeric(y)) {
    stop('y must be numeric (use 0/1 for binary response)')

  } else if (!is.numeric(tx) | length(setdiff(tx, 0:1)) > 0) {
    stop('tx must be vector of 0s and 1s')

  }

  #colnames(x) = paste('x', 1:ncol(x), sep = '')
  fit = list(x = x, pscore = pscore, ipcw = weight)

  z = tx * y / pscore - (1 - tx) * y / (1 - pscore)

  if (verbose) cat('fitting IPW treatment grf\n')

  data <- data.frame(y = z, x = x)
  colnames(data) <- c('y', colnames(x))

  if (meta_learner){
    fit$tau_fit <- regression_forest(as.matrix(data[,2:dim(data)[2]]),
                                     data$y,
                                     sample.weights = weight)
  }else{
    fit$tau_fit <- glm(y ~., family = "gaussian", data = data)
  }

  class(fit) = 'Fgrf'
  fit
}

predict.Fgrf = function(object, newx, meta_learner=TRUE, ...) {
  if (meta_learner){
    return(predict(object$tau_fit, newdata = as.matrix(newx))$predictions)
  }else{
    return(predict(object$tau_fit, newx))
  }
}
