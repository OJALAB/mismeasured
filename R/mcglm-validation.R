# ---------------------------------------------------------------------------
# Main study / validation study designs (Yi, Yan, Liao and Spiegelman, 2019,
# Section 4).
#
# A validation sample V observes the true category Z next to the proxy
# Z_hat. It is either internal (a subsample of the regression rows, so Y and
# x are also available for it) or external (separate units without Y). The
# misclassification nuisance eta is estimated from V by closed-form
# (multinomial ML) proportions; the regression parameter psi then solves
#
#   internal: sum_{i in M \ V} U_i(psi, eta) + sum_{i in V} S_i(psi) = 0,
#   external: sum_{i in M} U_i(psi, eta) = 0,
#
# where U is the method's corrected estimating function and S the score at
# the true Z. Inference stacks theta = (psi, eta) over all units:
#   Var(theta_hat) = A^{-1} B A^{-T},  A = sum_u w_u d g_u / d theta',
#   B = sum_u w_u g_u g_u',
# which covers both designs without an explicit n_V / n ratio.
#
# Nuisance parameterization (K categories, baseline category 0 implicit):
#   Pi[j, l], j = 1..K-1, l = 0..K-1    from V:   1{Z = l}(1{Z_hat = j} - Pi[j, l])
#   pi_source = "validation": pi_l, l = 1..K-1 from V:   1{Z = l} - pi_l
#   pi_source = "main":       p_j,  j = 1..K-1 from M:   1{Z_hat = j} - p_j,
#                             with pi = Pi^{-1} p (only Pi is transported).
# ---------------------------------------------------------------------------

# Methods whose misclassification nuisance is estimated from a validation
# sample (rather than supplied) when mcglm(validation = ) is used.
.mcglm_validated_methods <- c("sub")

#' Parse and check a validation-sample description
#'
#' @param validation \code{list(z, z_hat)} (external) or
#'   \code{list(z, index)} (internal).
#' @param z_hat Main-study proxy codes (0-based).
#' @param K Number of categories.
#' @return List with \code{z}, \code{proxy}, \code{index}, \code{n},
#'   \code{type}.
#' @keywords internal
.mcglm_parse_validation <- function(validation, z_hat, K) {
  if (!is.list(validation) || is.null(validation$z) ||
      any(!names(validation) %in% c("z", "z_hat", "index")))
    stop("validation must be a list with z and z_hat (external) or index (internal).",
         call. = FALSE)
  z <- validation$z
  nv <- length(z)
  check_codes <- function(v) {
    is.numeric(v) && is.null(dim(v)) && length(v) == nv &&
      all(is.finite(v)) && all(v == floor(v)) && all(v >= 0 & v < K)
  }
  if (nv < 2L || !check_codes(z))
    stop("validation$z must contain at least two codes in 0, ..., K-1.",
         call. = FALSE)
  index <- validation$index
  proxy <- validation$z_hat
  if (!is.null(index)) {
    if (!is.numeric(index) || !is.null(dim(index)) || length(index) != nv ||
        any(!is.finite(index)) || any(index != floor(index)) ||
        any(index < 1 | index > length(z_hat)) || anyDuplicated(index))
      stop("validation$index must contain distinct regression row numbers, one per z.",
           call. = FALSE)
    if (!is.null(proxy) &&
        (!check_codes(proxy) || any(proxy != z_hat[index])))
      stop("validation$z_hat does not match z_hat at validation$index.",
           call. = FALSE)
    proxy <- z_hat[index]
  }
  if (!check_codes(proxy))
    stop("validation$z_hat must contain one code in 0, ..., K-1 per z.",
         call. = FALSE)
  if (any(tabulate(z + 1L, nbins = K) == 0L))
    stop("Every true category must occur in validation.", call. = FALSE)
  list(z = as.integer(z), proxy = as.integer(proxy),
       index = if (is.null(index)) NULL else as.integer(index), n = nv,
       type = if (is.null(index)) "external" else "internal")
}

#' Empirical (Pi, pi_z) from a validation sample
#'
#' Column-conditional rates \eqn{\hat\Pi_{j\ell} = \#(\hat Z = j, Z = \ell)
#' / \#(Z = \ell)} and true-category proportions, in the format expected
#' by \code{mcglm} (scalar \code{pi_z} when \eqn{K = 2}).
#' @keywords internal
.mcglm_validation_probabilities <- function(vd, K) {
  B <- matrix(tabulate(vd$proxy + 1L + K * vd$z, nbins = K * K), K) / vd$n
  pi_v <- colSums(B)
  list(Pi = sweep(B, 2, pi_v, "/"),
       pi_z = if (K == 2L) pi_v[2L] else pi_v)
}

#' Misclassification nuisance estimated from a validation sample
#'
#' @param vd Parsed validation description.
#' @param z_hat Main-study proxy codes.
#' @param K Number of categories.
#' @param pi_source \code{"validation"} or \code{"main"}.
#' @param wt Main-study frequency weights (\code{NULL} for unit weights);
#'   internal validation rows inherit them.
#' @return List with \code{eta}, \code{map(eta)} returning
#'   \code{list(Pi, pi)}, per-unit estimating-function rows
#'   \code{rows_v} (validation units) and \code{rows_m} (main units), the
#'   sum-scale Jacobian \code{A} and validation-unit weights \code{w_v}.
#' @keywords internal
.mcglm_nuisance_setup <- function(vd, z_hat, K, pi_source, wt = NULL) {
  s  <- K - 1L
  n  <- length(z_hat)
  w_m <- if (is.null(wt)) rep(1, n) else wt
  w_v <- if (vd$type == "internal") w_m[vd$index] else rep(1, vd$n)

  Zv <- outer(vd$z, 0:s, "==") * 1          # n_V x K indicators of Z
  Hv <- outer(vd$proxy, 0:s, "==") * 1      # n_V x K indicators of Z_hat
  Hm <- outer(z_hat, 0:s, "==") * 1         # n x K

  zw <- colSums(w_v * Zv)                   # weighted count of each Z
  Pi_hat <- crossprod(Hv * w_v, Zv) / rep(zw, each = K)
  n_pi <- K * s
  if (pi_source == "validation") {
    prev <- colSums(w_v * Zv)[-1L] / sum(w_v)
  } else {
    prev <- colSums(w_m * Hm)[-1L] / sum(w_m)
  }
  eta <- c(as.numeric(Pi_hat[-1L, , drop = FALSE]), prev)
  q <- length(eta)

  map <- function(eta) {
    eta <- unname(eta)
    P <- matrix(eta[seq_len(n_pi)], s, K)
    Pi <- rbind(1 - colSums(P), P)
    prev <- eta[n_pi + seq_len(s)]
    if (pi_source == "validation") {
      pi <- c(1 - sum(prev), prev)
    } else {
      pi <- tryCatch(as.numeric(solve(Pi, c(1 - sum(prev), prev))),
                     error = function(e)
                       stop("pi_source = 'main': the estimated Pi is ",
                            "singular.", call. = FALSE))
    }
    list(Pi = Pi, pi = pi)
  }

  # Per-unit rows at eta_hat; columns ordered as eta.
  rows_v <- matrix(0, vd$n, q)
  for (l in seq_len(K)) {
    cols <- (l - 1L) * s + seq_len(s)
    rows_v[, cols] <- Zv[, l] *
      (Hv[, -1L, drop = FALSE] - rep(Pi_hat[-1L, l], each = vd$n))
  }
  rows_m <- matrix(0, n, q)
  prev_cols <- n_pi + seq_len(s)
  if (pi_source == "validation") {
    rows_v[, prev_cols] <- Zv[, -1L, drop = FALSE] - rep(prev, each = vd$n)
  } else {
    rows_m[, prev_cols] <- Hm[, -1L, drop = FALSE] - rep(prev, each = n)
  }

  # d(sum of rows)/d eta is diagonal: minus the weighted count entering each
  # proportion.
  A <- diag(c(rep(-zw, each = s),
              rep(if (pi_source == "validation") -sum(w_v) else -sum(w_m),
                  s)), q)

  pi_hat <- map(eta)$pi
  if (any(!is.finite(pi_hat)) || any(pi_hat <= 0 | pi_hat >= 1))
    stop("pi_source = 'main': Pi^{-1} applied to the main-study proxy ",
         "frequencies gives prevalences (",
         paste(round(pi_hat, 4), collapse = ", "),
         ") outside (0, 1); the validation Pi is inconsistent with the ",
         "main study. Use pi_source = 'validation'.", call. = FALSE)

  nms <- c(as.vector(outer(seq_len(s), 0:s,
                           function(j, l) sprintf("Pi[%d,%d]", j, l))),
           sprintf(if (pi_source == "validation") "pi[%d]" else "p_hat[%d]",
                   seq_len(s)))
  names(eta) <- nms
  list(eta = eta, map = map, rows_v = rows_v, rows_m = rows_m, A = A,
       w_v = w_v, w_m = w_m, pi_source = pi_source, names = nms)
}

#' Score at the true category for internal validation rows
#'
#' @return List with \code{U} (rows \eqn{\xi_i\{Y_i - \mu(\psi^\top\xi_i)\}})
#'   and \code{xi}, \code{mud} for the Jacobian.
#' @keywords internal
.mcglm_true_parts <- function(psi, y, z, x, K, fam) {
  xi  <- .mcglm_build_xi_hat(z, x, K)
  eta <- as.numeric(xi %*% psi)
  list(U = xi * (y - fam$linkinv(eta)), xi = xi, mud = fam$mu.eta(eta))
}

#' Sum-scale Jacobian of the true-score rows
#' @keywords internal
.mcglm_true_jacobian <- function(parts, w) {
  -crossprod(parts$xi * (w * parts$mud), parts$xi)
}

#' Stacked sandwich for (psi, eta) under a validation design
#'
#' @param psi Point estimate.
#' @param nuis Output of \code{.mcglm_nuisance_setup}.
#' @param vd Parsed validation description.
#' @param main_rows \code{function(psi, eta)} returning the n x p matrix of
#'   corrected estimating-function rows for every main-study unit.
#' @param main_jac \code{function(psi, eta, w)}: sum-scale Jacobian of
#'   \code{colSums(w * main_rows(psi, eta))} with respect to \code{psi}.
#' @param true_parts For internal validation, \code{.mcglm_true_parts} at
#'   \code{psi} for the validation rows; \code{NULL} otherwise.
#' @return List with \code{vcov} (full, (p + q) square), \code{psi}
#'   (p x p block) and \code{eta} (q x q block).
#' @keywords internal
.mcglm_validation_sandwich <- function(psi, nuis, vd, main_rows, main_jac,
                                       true_parts = NULL) {
  p <- length(psi)
  q <- length(nuis$eta)
  n <- nrow(nuis$rows_m)
  w_use <- nuis$w_m
  if (vd$type == "internal") w_use[vd$index] <- 0

  U <- main_rows(psi, nuis$eta)
  A_pp <- main_jac(psi, nuis$eta, w_use)
  if (!is.null(true_parts))
    A_pp <- A_pp + .mcglm_true_jacobian(true_parts, nuis$w_v)

  # d sum(w U) / d eta by central differences (q is at most K^2 - 1).
  f <- function(e) colSums(w_use * main_rows(psi, e))
  A_pe <- vapply(seq_len(q), function(k) {
    h <- 1e-6
    e1 <- e2 <- nuis$eta
    e1[k] <- e1[k] + h
    e2[k] <- e2[k] - h
    (f(e1) - f(e2)) / (2 * h)
  }, numeric(p))
  A_pe <- matrix(A_pe, p, q)

  A <- rbind(cbind(A_pp, A_pe), cbind(matrix(0, q, p), nuis$A))

  # Unit rows: main units first, then (external) validation units.
  G_m <- cbind(U * (w_use > 0), nuis$rows_m)
  if (vd$type == "internal") {
    G_m[vd$index, seq_len(p)] <- true_parts$U
    G_m[vd$index, p + seq_len(q)] <- G_m[vd$index, p + seq_len(q)] +
      nuis$rows_v
    G <- G_m
    w <- nuis$w_m
  } else {
    G <- rbind(G_m, cbind(matrix(0, vd$n, p), nuis$rows_v))
    w <- c(nuis$w_m, nuis$w_v)
  }
  B <- crossprod(G * w, G)

  A_inv <- solve(A)
  V <- A_inv %*% B %*% t(A_inv)
  list(vcov = V,
       psi = V[seq_len(p), seq_len(p), drop = FALSE],
       eta = V[p + seq_len(q), p + seq_len(q), drop = FALSE])
}
