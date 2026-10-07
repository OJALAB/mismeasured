# ---------------------------------------------------------------------------
# Shrinkage estimators of the misclassification matrix
#
# With many categories some true classes have few audited units and the
# column proportions C[, l] / n_l are noisy or contain zeros. The
# empirical-Bayes estimator shrinks each column towards a target T:
#   Pi[, l] = (C[, l] + c T[, l]) / (n_l + c),
# the posterior mean under a Dirichlet(c T[, l]) prior. The concentration c
# is fixed or chosen by maximising the Dirichlet-multinomial marginal
# likelihood. For a fixed c the shrinkage is O(1 / n_l), so first-order
# inference is that of the column proportions; the stacked estimating
# function C_jl + c T_jl - Pi_jl (n_l + c) is zero at the estimate.
# ---------------------------------------------------------------------------

#' Shrinkage target for the misclassification matrix
#'
#' @param C K x K (weighted) counts, rows proxy, columns true category.
#' @param target \code{"pooled"}, \code{"groups"}, \code{"loglinear"} or
#'   \code{"matrix"}.
#' @param groups Category-to-group map (length K) for \code{"groups"} and
#'   optionally \code{"loglinear"}.
#' @param target_matrix User target (K x K, columns summing to one).
#' @return K x K column-stochastic matrix.
#' @keywords internal
.mc_shrink_target <- function(C, target, groups = NULL, target_matrix = NULL) {
  K <- nrow(C)
  N <- sum(C)
  diag_rate <- sum(diag(C)) / N
  if (target == "matrix") {
    T <- as.matrix(target_matrix)
    if (!identical(dim(T), c(K, K)) || any(!is.finite(T)) || any(T < 0) ||
        any(abs(colSums(T) - 1) > 1e-8))
      stop("target_matrix must be a K x K matrix with columns summing to one.",
           call. = FALSE)
    return(T)
  }
  if (target == "pooled") {
    T <- matrix((1 - diag_rate) / (K - 1L), K, K)
    diag(T) <- diag_rate
    return(T)
  }
  if (target %in% c("groups", "loglinear") && !is.null(groups)) {
    if (length(groups) != K || anyNA(groups))
      stop("groups must assign each of the K categories to a group.",
           call. = FALSE)
    same <- outer(groups, groups, "==") & !diag(K)
  } else {
    same <- NULL
  }
  if (target == "groups") {
    if (is.null(same))
      stop("target = 'groups' needs control_mc(groups = ).", call. = FALSE)
    other <- !same & !diag(K)
    rate_same <- sum(C[same]) / N
    rate_other <- sum(C[other]) / N
    T <- matrix(0, K, K)
    for (l in seq_len(K)) {
      ns <- sum(same[, l])
      no <- sum(other[, l])
      T[l, l] <- diag_rate
      if (ns) T[same[, l], l] <- rate_same / ns
      if (no) T[other[, l], l] <- rate_other / no
      T[, l] <- T[, l] / sum(T[, l])
    }
    return(T)
  }
  # loglinear: quasi-independence with a diagonal (and same-group) effect
  cells <- expand.grid(j = factor(seq_len(K)), l = factor(seq_len(K)))
  cells$count <- as.numeric(C)
  cells$diag <- as.numeric(cells$j == cells$l)
  form <- count ~ j + l + diag
  if (!is.null(same)) {
    cells$same <- as.numeric(same)
    form <- count ~ j + l + diag + same
  }
  fit <- suppressWarnings(stats::glm(form, family = stats::poisson(),
                                     data = cells))
  T <- matrix(stats::fitted(fit), K, K)
  sweep(T, 2, colSums(T), "/")
}

#' Concentration maximising the Dirichlet-multinomial marginal likelihood
#' @keywords internal
.mc_shrink_concentration <- function(C, T) {
  n_l <- colSums(C)
  loglik <- function(logc) {
    cc <- exp(logc)
    a <- cc * T
    sum(lgamma(cc) - lgamma(n_l + cc)) +
      sum(lgamma(C + a) - lgamma(a), na.rm = TRUE)
  }
  exp(stats::optimize(loglik, c(log(1e-3), log(1e6)), maximum = TRUE)$maximum)
}

#' Shrinkage specification from validation counts
#'
#' @param C K x K counts normalised to the number of validation units.
#' @param control A \code{\link{control_mc}} object.
#' @return \code{NULL} for \code{estimator = "mle"}, otherwise a list with
#'   \code{c} (concentration), \code{T} (target), \code{target} and
#'   \code{estimator}.
#' @keywords internal
.mc_shrink_setup <- function(C, control) {
  if (control$estimator == "mle") return(NULL)
  K <- nrow(C)
  if (control$estimator == "dirichlet") {
    T <- matrix(1 / K, K, K)
    warning("estimator = 'dirichlet' shrinks towards the uniform matrix, ",
            "which is singular; keep alpha small relative to the counts or ",
            "use estimator = 'eb'.", call. = FALSE)
    return(list(c = K * control$alpha, T = T, target = "uniform",
                estimator = "dirichlet"))
  }
  T <- .mc_shrink_target(C, control$target, control$groups,
                         control$target_matrix)
  kap <- tryCatch(kappa(T, exact = TRUE), error = function(e) Inf)
  if (!is.finite(kap) || kap > control$kappa_max)
    stop("The shrinkage target (", control$target, ") is ill-conditioned ",
         "(condition number ", signif(kap, 3), "); shrinking towards it ",
         "would make Pi unidentifiable. Use another target.", call. = FALSE)
  cc <- if (identical(control$concentration, "ml"))
    .mc_shrink_concentration(C, T) else control$concentration
  list(c = cc, T = T, target = control$target, estimator = "eb")
}

#' Shrinkage penalty in the estimating function of the free Pi entries
#'
#' \eqn{c(T_{j\ell} - \Pi_{j\ell})} for \eqn{j \ge 1}, ordered as the free
#' entries (column-major over \eqn{\ell}); zero without shrinkage.
#' @keywords internal
.mc_shrink_penalty <- function(Pi, cc, Tm) {
  if (cc == 0) return(numeric(length(Pi) - ncol(Pi)))
  as.numeric(cc * (Tm[-1L, , drop = FALSE] - Pi[-1L, , drop = FALSE]))
}
