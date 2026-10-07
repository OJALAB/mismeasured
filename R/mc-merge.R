# ---------------------------------------------------------------------------
# Suggesting categories to merge when the classifier cannot separate them
#
# Merging categories a and b (in both the proxy and the true category)
# collapses the misclassification matrix exactly:
#   rows:    Pi'[c, l] = Pi[a, l] + Pi[b, l]
#   columns: Pi'[, c]  = (pi_a Pi'[, a] + pi_b Pi'[, b]) / (pi_a + pi_b),
#   prevalence pi'_c = pi_a + pi_b.
# Mutual confusion of a and b becomes correct classification, which raises
# the smallest singular value of Pi when the two columns were close.
# ---------------------------------------------------------------------------

#' Collapse a misclassification matrix over a partition of the categories
#'
#' @param Pi K x K column-stochastic matrix.
#' @param pi Prevalence (length K).
#' @param clusters List of integer index vectors partitioning \code{1:K}.
#' @return List with the collapsed \code{Pi} and \code{pi}.
#' @keywords internal
.mc_collapse <- function(Pi, pi, clusters) {
  G <- length(clusters)
  P <- matrix(0, G, G)
  pr <- vapply(clusters, function(h) sum(pi[h]), numeric(1))
  for (g in seq_len(G)) for (h in seq_len(G))
    P[g, h] <- sum(Pi[clusters[[g]], clusters[[h]], drop = FALSE] %*%
                     pi[clusters[[h]]]) / pr[h]
  list(Pi = P, pi = pr)
}

#' Likelihood-ratio test that two true classes share a proxy distribution
#'
#' G-test of homogeneity of the audit counts \code{a} and \code{b} (proxy
#' categories, possibly already merged, by two true classes).
#' @return p-value (1 when there is nothing to compare).
#' @keywords internal
.mc_homogeneity_p <- function(a, b) {
  M <- cbind(a, b)
  M <- M[rowSums(M) > 0, , drop = FALSE]
  if (nrow(M) < 2L || any(colSums(M) == 0)) return(1)
  E <- outer(rowSums(M), colSums(M)) / sum(M)
  G <- 2 * sum(ifelse(M > 0, M * log(M / E), 0))
  stats::pchisq(G, nrow(M) - 1L, lower.tail = FALSE)
}

#' Conditioning summaries of a misclassification matrix
#' @keywords internal
.mc_conditioning <- function(Pi) {
  sv <- svd(Pi)$d
  K <- nrow(Pi)
  akn <- if (K == 2L) abs(1 - Pi[2L, 1L] - Pi[1L, 2L]) else
    tryCatch(kappa(Pi[-1L, -1L, drop = FALSE] - Pi[-1L, 1L], exact = TRUE),
             error = function(e) Inf)
  list(sigma_min = min(sv),
       kappa = if (min(sv) > 0) max(sv) / min(sv) else Inf, akn = akn)
}

#' Suggest categories to merge
#'
#' When the classifier cannot tell some categories apart, the estimated
#' misclassification matrix is nearly singular and every correction is
#' unstable; shrinkage (\code{control_mc(estimator = "eb")}) cannot restore
#' information the classifier does not provide. \code{suggest_merge()}
#' searches greedily for the merges that improve the conditioning most.
#' Merging categories \eqn{a} and \eqn{b} -- in the proxy and in the true
#' category -- collapses \eqn{\Pi} exactly: rows are summed,
#' \eqn{\Pi'_{c\ell} = \Pi_{a\ell} + \Pi_{b\ell}}, and columns are averaged
#' with the prevalences, \eqn{\Pi'_{\cdot c} = (\pi_a\Pi'_{\cdot a} +
#' \pi_b\Pi'_{\cdot b}) / (\pi_a + \pi_b)}. At each step the pair whose
#' merge gives the largest smallest singular value (\code{criterion =
#' "sigma_min"}, recommended) or the smallest condition number
#' (\code{"kappa"}) is merged. Merging continues while the
#' \code{\link{diagnose_mc}} thresholds of the estimate fail \emph{or} the
#' chosen pair's proxy distributions are not significantly different
#' (likelihood-ratio test of homogeneity of their audit counts,
#' \eqn{p > \alpha}): sampling noise makes nearly identical columns of
#' \eqn{\hat\Pi} look separated, so the estimated conditioning alone can pass
#' while the effects of the two categories remain unidentified in practice.
#' It stops when neither holds, no pair is allowed, or two categories
#' remain.
#'
#' The function only suggests: merging replaces two regression effects by
#' one (a prevalence-weighted mixture, interpretable only if the categories
#' are similar), assumes that the merged categories are misclassified in a
#' similar way, and choosing a merge from the same audit is a data-driven
#' model choice that later standard errors do not reflect. Apply a merge
#' with \code{\link{merge_levels}} to the data and the validation sample and
#' refit.
#'
#' @param object An \code{"mc_estimate"} from \code{\link{estimate_mc}}.
#' @param groups Optional group of each category (in the estimate's level
#'   order); merges are then only allowed within a group (e.g. within an
#'   ISCO major group).
#' @param criterion \code{"sigma_min"} (default) or \code{"kappa"}.
#' @param max_merges Maximum number of merges (default: until the stopping
#'   rule holds or two categories remain).
#' @param alpha Significance level of the homogeneity test: a pair whose
#'   columns are not significantly different at this level is merged even
#'   when the diagnostics pass. Use \code{alpha = 1} to stop as soon as the
#'   diagnostics pass.
#' @return An object of class \code{"mc_merge"} with \code{steps} (a data
#'   frame: merged pair, new label, number of categories, smallest singular
#'   value, condition number, the AKN conditioning measure, smallest number
#'   of audited units per true class, the pair's mutual confusion share,
#'   the homogeneity-test p-value of the merged pair and whether the
#'   diagnostics pass), \code{mappings} (the label mapping after
#'   each step) and \code{mapping} (the last one).
#' @seealso \code{\link{merge_levels}}, \code{\link{diagnose_mc}},
#'   \code{vignette("sparse", "mismeasured")}
#' @examples
#' set.seed(3)
#' K <- 4
#' Pi <- matrix(c(0.90, 0.04, 0.03, 0.03,
#'                0.04, 0.47, 0.45, 0.04,    # categories b and c are
#'                0.04, 0.45, 0.47, 0.04,    # hardly distinguishable
#'                0.03, 0.03, 0.04, 0.90), K, K)
#' z <- sample(0:3, 2000, TRUE)
#' z_hat <- vapply(z, function(k) sample(0:3, 1, prob = Pi[, k + 1]), 0)
#' lev <- c("a", "b", "c", "d")
#' est <- estimate_mc(validation_sample(lev[z + 1], lev[z_hat + 1]),
#'                    control = control_mc(on_ill = "none"))
#' suggest_merge(est)
#' @export
suggest_merge <- function(object, groups = NULL,
                          criterion = c("sigma_min", "kappa"),
                          max_merges = NULL, alpha = 0.05) {
  if (!inherits(object, "mc_estimate"))
    stop("suggest_merge() needs an object from estimate_mc().", call. = FALSE)
  criterion <- match.arg(criterion)
  K <- object$K
  lev <- object$levels
  if (!is.null(groups) && (length(groups) != K || anyNA(groups)))
    stop("groups must give a group for each of the ", K, " categories.",
         call. = FALSE)
  ctrl <- object$control
  Pi <- unname(object$Pi)
  pi <- unname(object$pi)
  vb <- object$validation
  C <- matrix(tabulate(vb$proxy + 1L + K * vb$z, nbins = K * K), K, K)
  n_class <- colSums(C)
  if (is.null(max_merges)) max_merges <- K - 2L

  clusters <- as.list(seq_len(K))
  label_of <- function(h) paste(lev[h], collapse = "+")
  mapping_of <- function(cl) {
    out <- character(K)
    for (h in cl) out[h] <- label_of(h)
    stats::setNames(out, lev)
  }
  summarise <- function(cl) {
    col <- .mc_collapse(Pi, pi, cl)
    cond <- .mc_conditioning(col$Pi)
    n_min <- min(vapply(cl, function(h) sum(n_class[h]), numeric(1)))
    akn_ok <- if (length(cl) == 2L) cond$akn >= ctrl$sigma_min else
      cond$akn <= ctrl$kappa_max
    c(cond, list(n_min = n_min,
                 pass = cond$sigma_min >= ctrl$sigma_min &&
                   cond$kappa <= ctrl$kappa_max && akn_ok &&
                   n_min >= ctrl$min_class_n))
  }
  row_of <- function(step, merged, label, cl, share, p) {
    sm <- summarise(cl)
    data.frame(step = step, merged = merged, new_label = label,
               K = length(cl), sigma_min = sm$sigma_min, kappa = sm$kappa,
               akn = sm$akn, min_class_n = sm$n_min,
               confusion_share = share, p_homogeneity = p, pass = sm$pass,
               stringsAsFactors = FALSE)
  }
  # Proxy counts of a true cluster, with proxy rows merged like the clusters.
  counts_of <- function(h, cl)
    vapply(cl, function(g) sum(C[g, h]), numeric(1))

  steps <- row_of(0L, NA_character_, NA_character_, clusters, NA_real_,
                  NA_real_)
  mappings <- list(mapping_of(clusters))
  step <- 0L
  while (length(clusters) > 2L && step < max_merges) {
    best <- NULL
    best_value <- if (criterion == "sigma_min") -Inf else Inf
    G <- length(clusters)
    for (a in seq_len(G - 1L)) for (b in (a + 1L):G) {
      if (!is.null(groups) &&
          groups[clusters[[a]][1L]] != groups[clusters[[b]][1L]]) next
      cand <- clusters[-c(a, b)]
      cand <- c(cand, list(sort(c(clusters[[a]], clusters[[b]]))))
      cand <- cand[order(vapply(cand, min, numeric(1)))]
      cond <- .mc_conditioning(.mc_collapse(Pi, pi, cand)$Pi)
      value <- if (criterion == "sigma_min") cond$sigma_min else cond$kappa
      better <- if (criterion == "sigma_min") value > best_value else
        value < best_value
      if (better) {
        best_value <- value
        best <- list(a = a, b = b, cand = cand)
      }
    }
    if (is.null(best)) break
    ia <- clusters[[best$a]]
    ib <- clusters[[best$b]]
    p_hom <- .mc_homogeneity_p(counts_of(ia, clusters),
                               counts_of(ib, clusters))
    if (steps$pass[nrow(steps)] && p_hom <= alpha) break
    step <- step + 1L
    share <- (sum(C[ia, ib]) + sum(C[ib, ia])) /
      (sum(n_class[ia]) + sum(n_class[ib]))
    clusters <- best$cand
    steps <- rbind(steps, row_of(step, paste(label_of(ia), label_of(ib),
                                             sep = " & "),
                                 label_of(sort(c(ia, ib))), clusters, share,
                                 p_hom))
    mappings[[length(mappings) + 1L]] <- mapping_of(clusters)
  }
  structure(list(steps = steps, mappings = mappings,
                 mapping = mappings[[length(mappings)]],
                 criterion = criterion, groups = groups),
            class = "mc_merge")
}

#' @export
print.mc_merge <- function(x, digits = 3L, ...) {
  cat("Suggested merges (criterion: ", x$criterion,
      if (!is.null(x$groups)) "; within groups", ")\n", sep = "")
  s <- x$steps
  s$sigma_min <- signif(s$sigma_min, digits)
  s$kappa <- signif(s$kappa, digits)
  s$akn <- signif(s$akn, digits)
  s$confusion_share <- round(s$confusion_share, digits)
  s$p_homogeneity <- signif(s$p_homogeneity, digits)
  print(s, row.names = FALSE)
  if (nrow(x$steps) > 1L) {
    cat("\nMapping after the last step:\n")
    changed <- x$mapping[names(x$mapping) != x$mapping]
    print(changed)
  }
  if (!x$steps$pass[nrow(x$steps)])
    cat("\nThe diagnostics still flag a problem after the last merge.\n")
  cat("\nSuggestions only: a merge replaces two effects by one and is a",
      "data-driven model choice.\nApply it with merge_levels() and refit.\n")
  invisible(x)
}

#' Recode categories according to a merge
#'
#' Applies a mapping of category labels -- from \code{\link{suggest_merge}}
#' or a named character vector \code{c(old = "new")} -- to a factor or
#' character vector (e.g. the \code{mc()} variable in the data) or to a
#' \code{\link{validation_sample}} (its true and proxy categories, and the
#' stratum sizes of a \code{strata = "z_hat"} design). Labels not in the
#' mapping are kept. The merged levels keep the order of their first
#' original level, so the baseline category stays first.
#'
#' @param x A factor, a character vector or an \code{"mc_validation"}
#'   object.
#' @param mapping An \code{"mc_merge"} object or a named character vector.
#' @param step For an \code{"mc_merge"} object, the step whose mapping is
#'   used (default: the last).
#' @return The recoded factor, or the recoded validation sample.
#' @examples
#' merge_levels(factor(c("a", "b", "c", "b")), c(b = "b+c", c = "b+c"))
#' @export
merge_levels <- function(x, mapping, step = NULL) {
  if (inherits(mapping, "mc_merge")) {
    mapping <- if (is.null(step)) mapping$mapping else
      mapping$mappings[[step + 1L]]
  }
  if (!is.character(mapping) || is.null(names(mapping)))
    stop("mapping must be an mc_merge object or a named character vector.",
         call. = FALSE)
  recode <- function(v) {
    v <- as.character(v)
    hit <- v %in% names(mapping)
    v[hit] <- mapping[v[hit]]
    v
  }
  if (inherits(x, "mc_validation")) {
    out <- x
    out$z <- recode(x$z)
    if (!is.null(x$z_hat)) out$z_hat <- recode(x$z_hat)
    if (!is.null(x$label_order))
      out$label_order <- unique(recode(x$label_order))
    out$numeric_codes <- x$numeric_codes &&
      all(grepl("^[0-9]+$", c(out$z, out$z_hat)))
    if (x$strata_type == "z_hat" && !is.null(x$N_strata)) {
      Ns <- tapply(x$N_strata, recode(names(x$N_strata)), sum)
      out$N_strata <- stats::setNames(as.numeric(Ns), names(Ns))
    }
    return(out)
  }
  old <- if (is.factor(x)) levels(x) else sort(unique(as.character(x)))
  factor(recode(x), levels = unique(recode(old)))
}
