# ---------------------------------------------------------------------------
# Internal: GLM family normalization
# ---------------------------------------------------------------------------

#' Normalize a GLM family specification
#'
#' @param family A \code{\link[stats]{family}} object, family function, or
#'   character string naming a family function.
#' @param env Environment in which to resolve a character family name.
#' @return A standard \code{family} object.
#' @keywords internal
.normalize_family <- function(family, env = parent.frame()) {
  if (is.character(family)) {
    if (length(family) != 1L)
      stop("'family' must be a family object, family function, or one name.",
           call. = FALSE)
    family <- get(family, mode = "function", envir = env)()
  }
  if (is.function(family)) family <- family()
  if (!inherits(family, "family"))
    stop("'family' not recognized", call. = FALSE)

  .check_canonical_link(family)
  family
}
