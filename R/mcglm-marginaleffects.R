# ---------------------------------------------------------------------------
# marginaleffects integration for mcglm fits.
#
# marginaleffects gates models by class; .onLoad() (zzz.R) appends "mcglm"
# to the 'marginaleffects_model_classes' option, and the four methods
# below are registered lazily on the marginaleffects generics
# (marginaleffects in Suggests). They are internal hooks -- users just
# call the ordinary marginaleffects functions on a formula-interface
# fit, with no extra arguments:
#   avg_slopes(fit)
#   avg_comparisons(fit, variables = "z")
# (the original data are recovered from the stored call; pass newdata =
# to override). All quantities use the fit's default method (the last
# one fit, e.g. "cs"); the delta-method standard errors propagate
# vcov(fit) of that method.
# ---------------------------------------------------------------------------

#' @exportS3Method marginaleffects::get_coef
get_coef.mcglm <- function(model, ...) {
  stats::coef(model)
}

#' @exportS3Method marginaleffects::set_coef
set_coef.mcglm <- function(model, coefs, ...) {
  m <- .mcglm_default_method(model)
  model$coefficients[[m]][] <- coefs
  model
}

#' @exportS3Method marginaleffects::get_vcov
get_vcov.mcglm <- function(model, vcov = NULL, ...) {
  if (is.matrix(vcov)) return(vcov)
  if (is.function(vcov)) return(vcov(model))
  stats::vcov(model)
}

#' @exportS3Method marginaleffects::get_predict
get_predict.mcglm <- function(model, newdata, ...) {
  if (is.null(model$formula))
    stop("marginaleffects support requires an mcglm fit from the formula ",
         "interface (predictions on modified data need the mc() formula).")
  data.frame(
    rowid    = seq_len(nrow(newdata)),
    estimate = stats::predict(model, newdata = as.data.frame(newdata),
                              type = "response")
  )
}
