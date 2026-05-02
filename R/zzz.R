#' @useDynLib mismeasured, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Suppress R CMD check NOTE for NSE variables used in glm(weights = .wt)
utils::globalVariables(".wt")
