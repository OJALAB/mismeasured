# ---------------------------------------------------------------------------
# Validation-sample class
#
# A validation (audit) sample observes the true category Z next to the
# proxy Z_hat. It is internal (a subsample of the main-study rows, given by
# `index`) or external (separate units, proxies given by `z_hat`). Design
# information -- inverse inclusion probabilities, strata, frame sizes and
# the choice between Hajek and Horvitz-Thompson estimators -- travels with
# the sample so that estimate_mc() and mcglm() can use it.
#
# Categories are matched by label, never by position: factor or character
# values are compared with the model's category labels (the levels of the
# mc() variable), and numeric values are compared as labels too, so that a
# numeric audit matches an integer-coded model ("0", "1", ...) or a factor
# whose levels are the same numbers.
# ---------------------------------------------------------------------------

#' Validation sample for misclassification probabilities
#'
#' Describes a sample in which the true category \eqn{Z} is observed next
#' to its proxy \eqn{\hat Z}, together with its sampling design. The object
#' is passed to \code{\link{estimate_mc}} or directly to
#' \code{\link{mcglm}(validation = )}.
#'
#' @section Categories:
#' \code{z} and \code{z_hat} hold category labels: a factor, a character
#' vector or numbers. They are matched \emph{by label} to the categories of
#' the model -- the levels of the \code{mc()} variable in
#' \code{\link{mcglm}}, or \code{levels(z_hat)} / \code{levels} in
#' \code{\link{estimate_mc}} -- so the order of factor levels in the audit
#' does not matter. Numbers are compared as labels as well: they fit an
#' integer-coded model (categories \code{"0"}, \code{"1"}, ...) or a factor
#' whose levels are those numbers, and are rejected for a model with other
#' labels.
#'
#' @param z True categories of the validation units (factor, character or
#'   numeric labels).
#' @param z_hat Proxy categories for the same units (external validation).
#'   For internal validation they are taken from the main study at
#'   \code{index}; if supplied too, they must agree.
#' @param index Row numbers of the validated units in the main study
#'   (internal validation). \code{NULL} for an external sample.
#' @param weights Optional design weights \eqn{d_i} (inverse inclusion
#'   probabilities) of the validation units. For internal validation these
#'   are inclusion probabilities within the main study, e.g. for an audit
#'   subsample of a nonprobability sample. When \code{NULL}, weights follow
#'   from the design: \eqn{N_h / n_h} within strata, or a constant for a
#'   simple random sample.
#' @param strata Optional stratum labels of the validation units, or the
#'   string \code{"z_hat"} for an audit stratified by the predicted
#'   category. Strata define the default weights and enter the
#'   (stratum-centred) variance of an external sample.
#' @param estimator \code{"hajek"} (default; ratio estimators) or
#'   \code{"ht"} (Horvitz--Thompson; divides weighted totals by the known
#'   frame size \code{N} or stratum sizes \code{N_strata}). Column ratios of
#'   \eqn{\Pi} are identical under both; they differ for the prevalence.
#' @param N Size of the frame the sample was drawn from. Defaults to the
#'   (weighted) number of main-study rows for internal validation; required
#'   for \code{estimator = "ht"} with an external sample.
#' @param N_strata Named numeric vector of stratum sizes in the frame; for
#'   \code{strata = "z_hat"} the names are the category labels. Computed
#'   from the main-study proxies for internal validation with
#'   \code{strata = "z_hat"}; otherwise required when \code{strata} is
#'   given without \code{weights}.
#'
#' @return An object of class \code{"mc_validation"}.
#' @seealso \code{\link{estimate_mc}}, \code{\link{mcglm}}
#' @examples
#' set.seed(1)
#' z     <- rbinom(200, 1, 0.4)
#' z_hat <- ifelse(z == 1, rbinom(200, 1, 0.85), rbinom(200, 1, 0.1))
#' val <- validation_sample(z, z_hat)
#' val
#'
#' # labels work directly; their order does not matter
#' occ <- c("clerk", "manager")
#' validation_sample(factor(occ[z + 1]), occ[z_hat + 1])
#' @export
validation_sample <- function(z, z_hat = NULL, index = NULL, weights = NULL,
                              strata = NULL, estimator = c("hajek", "ht"),
                              N = NULL, N_strata = NULL) {
  estimator <- match.arg(estimator)
  nv <- length(z)
  is_labels <- function(v) {
    is.null(dim(v)) && length(v) == nv && !anyNA(v) &&
      (is.factor(v) || is.character(v) ||
         (is.numeric(v) && all(is.finite(v)) && all(v == floor(v))))
  }
  if (nv < 2L || !is_labels(z))
    stop("validation$z must contain at least two category labels or codes ",
         "(factor, character or whole numbers) without missing values.",
         call. = FALSE)
  if (is.null(z_hat) && is.null(index))
    stop("validation needs z_hat (external) or index (internal).",
         call. = FALSE)
  if (!is.null(z_hat) && !is_labels(z_hat))
    stop("validation$z_hat must contain one category label or code per z.",
         call. = FALSE)
  if (!is.null(index) &&
      (!is.numeric(index) || !is.null(dim(index)) || length(index) != nv ||
       any(!is.finite(index)) || any(index != floor(index)) ||
       any(index < 1) || anyDuplicated(index)))
    stop("validation$index must contain distinct regression row numbers, one per z.",
         call. = FALSE)
  if (!is.null(weights) &&
      (!is.numeric(weights) || length(weights) != nv ||
       any(!is.finite(weights)) || any(weights <= 0)))
    stop("validation weights must be ", nv, " finite positive design ",
         "weights (inverse inclusion probabilities).", call. = FALSE)
  strata_type <- "none"
  if (!is.null(strata)) {
    if (identical(strata, "z_hat")) {
      strata_type <- "z_hat"
    } else {
      if (length(strata) != nv || anyNA(strata))
        stop("strata must be \"z_hat\" or ", nv, " non-missing labels.",
             call. = FALSE)
      strata_type <- "user"
    }
  }
  if (!is.null(N) &&
      (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0))
    stop("N must be one positive frame size.", call. = FALSE)
  if (!is.null(N_strata) &&
      (!is.numeric(N_strata) || is.null(names(N_strata)) ||
       any(!is.finite(N_strata)) || any(N_strata <= 0)))
    stop("N_strata must be a named vector of positive stratum sizes.",
         call. = FALSE)
  if (estimator == "ht" && is.null(index) && is.null(N) &&
      is.null(N_strata))
    stop("estimator = 'ht' with an external validation sample needs the ",
         "frame size N (or N_strata).", call. = FALSE)

  as_label <- function(v) if (is.numeric(v)) as.character(as.integer(v)) else
    as.character(v)
  label_order <- unique(c(if (is.factor(z)) levels(z),
                          if (is.factor(z_hat)) levels(z_hat)))
  structure(list(z = as_label(z),
                 z_hat = if (is.null(z_hat)) NULL else as_label(z_hat),
                 numeric_codes = is.numeric(z) &&
                   (is.null(z_hat) || is.numeric(z_hat)),
                 label_order = if (length(label_order)) label_order else NULL,
                 index = if (is.null(index)) NULL else as.integer(index),
                 weights = if (is.null(weights)) NULL else as.numeric(weights),
                 strata = if (strata_type == "user") strata else NULL,
                 strata_type = strata_type, estimator = estimator,
                 N = N, N_strata = N_strata, n = nv,
                 type = if (is.null(index)) "external" else "internal"),
            class = "mc_validation")
}

#' Coerce to a validation sample
#'
#' Converts the list form accepted by \code{mcglm(validation = )} --
#' \code{list(z, z_hat)} (external) or \code{list(z, index)} (internal) --
#' into a \code{\link{validation_sample}}.
#' @param x A list or an \code{"mc_validation"} object.
#' @return An object of class \code{"mc_validation"}.
#' @export
as_validation_sample <- function(x) {
  if (inherits(x, "mc_validation")) return(x)
  if (!is.list(x) || is.null(x$z) ||
      any(!names(x) %in% c("z", "z_hat", "index")))
    stop("validation must be a list with z and z_hat (external) or index (internal).",
         call. = FALSE)
  validation_sample(x$z, z_hat = x$z_hat, index = x$index)
}

#' @export
print.mc_validation <- function(x, ...) {
  design <- c(if (!is.null(x$weights)) "design weights",
              switch(x$strata_type, z_hat = "stratified by z_hat",
                     user = "stratified", none = NULL))
  cat(sprintf("Validation sample (%s): n = %d, %s, %s estimator\n",
              x$type, x$n,
              if (length(design)) paste(design, collapse = ", ") else
                "simple random sample",
              if (x$estimator == "ht") "Horvitz-Thompson" else "Hajek"))
  if (!is.null(x$z_hat)) {
    print(table(z_hat = x$z_hat, z = x$z))
  } else {
    cat("True categories:\n")
    print(table(z = x$z))
    cat("Proxies are taken from main-study rows at `index`.\n")
  }
  invisible(x)
}

#' Category labels shared by a validation sample and a main study
#'
#' Determines the model's category labels and converts the main-study
#' proxies to 0-based codes.
#'
#' @param val An \code{"mc_validation"} object.
#' @param z_hat Main-study proxies (factor, character, integer codes) or
#'   \code{NULL}.
#' @param K Optional number of categories.
#' @param levels Optional category labels (fix order and K).
#' @return List with \code{levels} and \code{z_hat} (0-based codes or
#'   \code{NULL}).
#' @keywords internal
.mc_category_levels <- function(val, z_hat = NULL, K = NULL, levels = NULL) {
  if (!is.null(levels)) {
    levels <- as.character(levels)
    if (anyDuplicated(levels) || anyNA(levels))
      stop("levels must be distinct, non-missing category labels.",
           call. = FALSE)
  } else if (is.factor(z_hat)) {
    levels <- base::levels(z_hat)
  } else if (is.character(z_hat)) {
    levels <- base::levels(factor(z_hat))
  } else if (!is.null(z_hat) || val$numeric_codes) {
    codes <- c(if (!is.null(z_hat)) {
      .mcglm_check_z_hat(z_hat)
      as.integer(z_hat)
    }, if (val$numeric_codes) as.integer(c(val$z, val$z_hat)))
    if (any(codes < 0L))
      stop("Integer category codes must be 0, ..., K-1.", call. = FALSE)
    Kc <- max(2L, max(codes) + 1L, if (is.null(K)) 0L else K)
    levels <- as.character(seq_len(Kc) - 1L)
  } else {
    levels <- unique(c(val$label_order,
                       sort(unique(c(val$z, val$z_hat)))))
  }
  if (!is.null(K) && K != length(levels))
    stop("K = ", K, " does not match the ", length(levels),
         " category labels (", paste(levels, collapse = ", "), ").",
         call. = FALSE)
  codes <- NULL
  if (!is.null(z_hat)) {
    codes <- match(as.character(z_hat), levels) - 1L
    if (anyNA(codes))
      stop("z_hat has value(s) ",
           paste0("'", unique(as.character(z_hat)[is.na(codes)]), "'",
                  collapse = ", "),
           " that are not among the categories (",
           paste(levels, collapse = ", "), ").", call. = FALSE)
  }
  list(levels = levels, z_hat = codes)
}

#' Bind a validation sample to the main study
#'
#' Maps the audit's labels to the model's category codes, checks them
#' against the main-study proxies, resolves internal proxies, and computes
#' the design weights, frame sizes and strata used by the estimating
#' functions.
#'
#' @param val An \code{"mc_validation"} object.
#' @param z_hat Main-study proxies as 0-based codes (required for internal
#'   validation and for main-study-based quantities; may be \code{NULL}
#'   otherwise).
#' @param K Number of categories.
#' @param wt Main-study frequency weights (\code{NULL} for unit weights).
#' @param levels Category labels of the model (code k has label
#'   \code{levels[k + 1]}); defaults to \code{"0", ..., "K-1"}.
#' @return List with \code{type}, \code{z}, \code{proxy} (0-based codes),
#'   \code{index}, \code{n}, \code{d} (design weights), \code{w_v}
#'   (frequency weights of validation rows), \code{strata} (factor or
#'   \code{NULL}), \code{N}, \code{N_h} (frame size per validation unit's
#'   stratum), \code{estimator}, \code{user_weights} and \code{levels}.
#' @keywords internal
.mc_bind_validation <- function(val, z_hat, K, wt = NULL,
                                levels = as.character(seq_len(K) - 1L)) {
  nv <- val$n
  to_code <- function(lab, what) {
    code <- match(lab, levels) - 1L
    if (anyNA(code)) {
      numeric_levels <- identical(levels, as.character(seq_along(levels) - 1L))
      hint <- if (val$numeric_codes && !numeric_levels)
        paste0(" Numbers are matched to the category labels; give the ",
               "audit as labels, e.g. factor(z, levels = levels(data$v)) ",
               "for the mc() variable v.") else ""
      stop(what, " has label(s) ",
           paste0("'", unique(lab[is.na(code)]), "'", collapse = ", "),
           " that are not categories of the model (",
           paste(levels, collapse = ", "), ").", hint, call. = FALSE)
    }
    code
  }
  z <- to_code(val$z, "validation$z")
  proxy <- if (is.null(val$z_hat)) NULL else
    to_code(val$z_hat, "validation$z_hat")
  index <- val$index
  if (!is.null(index)) {
    if (is.null(z_hat))
      stop("Internal validation needs the main-study proxies z_hat.",
           call. = FALSE)
    if (any(index > length(z_hat)))
      stop("validation$index must contain distinct regression row numbers, one per z.",
           call. = FALSE)
    if (!is.null(proxy) && any(proxy != z_hat[index]))
      stop("validation$z_hat does not match z_hat at validation$index.",
           call. = FALSE)
    proxy <- as.integer(z_hat[index])
  }
  missing_true <- which(tabulate(z + 1L, nbins = K) == 0L)
  if (length(missing_true))
    stop("Every true category must occur in validation (none for ",
         paste0("'", levels[missing_true], "'", collapse = ", "), ").",
         call. = FALSE)

  internal <- !is.null(index)
  wt_m <- if (is.null(wt)) rep(1, length(z_hat)) else wt
  if (internal && !is.null(val$weights) && any(wt_m != 1))
    stop("Design weights for an internal validation sample cannot be ",
         "combined with main-study frequency weights.", call. = FALSE)
  w_v <- if (internal) wt_m[index] else rep(1, nv)

  strata <- switch(val$strata_type,
                   z_hat = factor(levels[proxy + 1L], levels = levels),
                   user = factor(val$strata),
                   none = NULL)
  N <- val$N
  if (is.null(N) && internal) N <- sum(wt_m)

  # Frame size of each validation unit's stratum (for default weights and
  # for HT). For z_hat strata of an internal sample it comes from the main
  # study's proxy counts.
  N_h <- NULL
  if (!is.null(strata)) {
    lev <- base::levels(droplevels(strata))
    Nvec <- val$N_strata
    if (is.null(Nvec) && val$strata_type == "z_hat" && internal) {
      Nvec <- stats::setNames(
        as.numeric(tapply(wt_m, factor(z_hat, levels = seq_len(K) - 1L),
                          sum)),
        levels)
    }
    if (!is.null(Nvec)) {
      if (!all(lev %in% names(Nvec)))
        stop("N_strata must name every stratum (",
             paste(lev, collapse = ", "), ").", call. = FALSE)
      N_h <- as.numeric(Nvec[as.character(strata)])
    } else if (is.null(val$weights)) {
      stop("Stratified validation without weights needs N_strata.",
           call. = FALSE)
    }
  }

  user_weights <- !is.null(val$weights)
  d <- if (user_weights) {
    val$weights
  } else if (!is.null(N_h)) {
    n_h <- as.numeric(table(strata)[as.character(strata)])
    N_h / n_h
  } else if (!is.null(N)) {
    rep(N / sum(w_v), nv)
  } else {
    rep(1, nv)
  }

  list(type = if (internal) "internal" else "external",
       z = z, proxy = proxy, index = index, n = nv, d = d, w_v = w_v,
       strata = strata, N = N, N_h = N_h, estimator = val$estimator,
       user_weights = user_weights, levels = levels)
}
