#' @useDynLib mismeasured, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Suppress R CMD check NOTE for NSE variables used in glm(weights = .wt)
utils::globalVariables(".wt")

.onLoad <- function(libname, pkgname) {
  # marginaleffects gates models by class; opt mcglm in (the actual
  # methods live in mcglm-marginaleffects.R). Append, don't clobber.
  cls <- getOption("marginaleffects_model_classes", default = NULL)
  if (!"mcglm" %in% cls)
    options(marginaleffects_model_classes = c(cls, "mcglm"))
  invisible()
}
