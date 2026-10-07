# ---------------------------------------------------------------------------
# Diagnostics for estimated misclassification matrices
# ---------------------------------------------------------------------------

#' Diagnose an estimated misclassification matrix
#'
#' Reports quantities that signal an unreliable correction: the smallest
#' singular value and condition number of \eqn{\hat\Pi}, the number of
#' validation units per true category, empty cells, the conditioning of the
#' matrix inverted by the AKN corrected score (\code{method = "cs_akn"}),
#' and whether the estimated prevalence lies on the boundary.
#'
#' @param object An \code{"mc_estimate"} from \code{\link{estimate_mc}}.
#' @param kappa_max,sigma_min,min_class_n Thresholds; default to those in
#'   the estimate's \code{\link{control_mc}}.
#' @return An object of class \code{"mc_diagnostics"}; its
#'   \code{problems} element lists the thresholds that are violated.
#' @seealso \code{\link{estimate_mc}}
#' @export
diagnose_mc <- function(object,
                        kappa_max = object$control$kappa_max,
                        sigma_min = object$control$sigma_min,
                        min_class_n = object$control$min_class_n) {
  if (!inherits(object, "mc_estimate"))
    stop("diagnose_mc() needs an object from estimate_mc().", call. = FALSE)
  K  <- object$K
  Pi <- object$Pi
  n_class <- tabulate(object$validation$z + 1L, nbins = K)
  sv <- svd(Pi)$d
  smin <- min(sv)
  kap <- if (smin > 0) max(sv) / smin else Inf
  if (K == 2L) {
    akn <- c(`|1 - p01 - p10|` = abs(1 - Pi[2, 1] - Pi[1, 2]))
  } else {
    Q <- Pi[-1L, -1L, drop = FALSE] - Pi[-1L, 1L]
    akn <- c(`kappa(Q)` = tryCatch(kappa(Q, exact = TRUE),
                                   error = function(e) Inf))
  }
  boundary <- isTRUE(object$prevalence$boundary)

  problems <- character()
  if (smin < sigma_min)
    problems <- c(problems, sprintf("smallest singular value of Pi %.3g < %g",
                                    smin, sigma_min))
  if (kap > kappa_max)
    problems <- c(problems, sprintf("condition number of Pi %.3g > %g",
                                    kap, kappa_max))
  if (any(n_class < min_class_n))
    problems <- c(problems, sprintf(
      "true categor%s %s with fewer than %d validation units",
      if (sum(n_class < min_class_n) > 1L) "ies" else "y",
      paste0("'", (if (is.null(object$levels)) seq_len(K) - 1L else
        object$levels)[n_class < min_class_n], "'", collapse = ", "),
      as.integer(min_class_n)))
  if (boundary)
    problems <- c(problems, "estimated prevalence on the boundary (variance unreliable)")

  lev <- if (is.null(object$levels)) as.character(seq_len(K) - 1L) else
    object$levels
  structure(list(n_class = stats::setNames(n_class, lev),
                 zero_cells = sum(Pi == 0), sigma_min = smin, kappa = kap,
                 akn = akn, boundary = boundary, problems = problems),
            class = "mc_diagnostics")
}

#' @export
print.mc_diagnostics <- function(x, ...) {
  cat("Misclassification diagnostics\n")
  cat(sprintf("  smallest singular value of Pi: %.4g\n", x$sigma_min))
  cat(sprintf("  condition number of Pi:        %.4g\n", x$kappa))
  cat(sprintf("  %s: %.4g\n", names(x$akn), x$akn))
  cat("  validation units per true category:",
      paste(names(x$n_class), x$n_class, sep = "=", collapse = ", "), "\n")
  cat("  empty cells of Pi:", x$zero_cells, "\n")
  if (length(x$problems)) {
    cat("  Problems:\n")
    for (p in x$problems) cat("   -", p, "\n")
  } else {
    cat("  No problems detected.\n")
  }
  invisible(x)
}

#' Act on diagnostic problems
#' @keywords internal
.mc_on_ill <- function(diag, on_ill) {
  if (!length(diag$problems) || on_ill == "none") return(invisible(NULL))
  if (on_ill == "regularize") on_ill <- "warn"
  msg <- paste0("Estimated misclassification matrix: ",
                paste(diag$problems, collapse = "; "), ".")
  if (on_ill == "error") stop(msg, call. = FALSE)
  warning(msg, call. = FALSE)
  invisible(NULL)
}
