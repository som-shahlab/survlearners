library(testthat)
library(survlearners)

# This treats warnings as errors
options(warn = 2)
test_check("survlearners")
