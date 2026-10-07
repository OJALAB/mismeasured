# ---------------------------------------------------------------------------
# Estimated misclassification probabilities
#
# estimate_mc() turns a validation sample into estimates of the
# column-conditional misclassification matrix Pi[j, l] = P(Z_hat = j | Z = l),
# the latent prevalence pi and the predictive matrix
# W[j, l] = P(Z = l | Z_hat = j), with a sandwich covariance from their
# (design-weighted) estimating functions. The same estimating functions are
# stacked with a regression estimator's in mcglm() to propagate their
# uncertainty.
#
# Nuisance vector eta (K categories, baseline 0 implicit):
#   Pi[j, l], j = 1..K-1, l = 0..K-1:
#       validation unit: d 1{Z = l} (1{Z_hat = j} - Pi[j, l])
#   prevalence = "validation", Hajek: pi_l, l >= 1
#       validation unit: d (1{Z = l} - pi_l)
#   prevalence = "validation", HT: pi_l, l >= 1
#       validation unit: d 1{Z = l}; total minus N pi_l
#       (internal: -pi_l on every main unit; external: -(N / n_V) pi_l on
#       every validation unit)
#   prevalence = "em": pi_l, l >= 1
#       main unit: P(Z = l | Z_hat_i; Pi, pi) - pi_l (EM fixed point)
#   prevalence = "inverse": p_j, j >= 1 with pi = Pi^{-1} p
#       main unit: 1{Z_hat_i = j} - p_j
# d are design weights; main and internal validation rows are additionally
# multiplied by main-study frequency weights in sums and meats.
# ---------------------------------------------------------------------------

#' Control settings for estimated misclassification probabilities
#'
#' @param estimator Estimator of \eqn{\Pi}: \code{"mle"} (design-weighted
#'   column proportions), \code{"eb"} (empirical-Bayes shrinkage of each
#'   column towards a target, \eqn{\hat\Pi_{\cdot\ell} = (C_{\cdot\ell} +
#'   c\,T_{\cdot\ell}) / (n_\ell + c)}) or \code{"dirichlet"} (add
#'   \code{alpha} to every cell; shrinks towards the singular uniform
#'   matrix). Shrinkage helps with many categories and few audited units
#'   per true class; it changes the estimate by \eqn{O(1/n_\ell)}, so the
#'   delta-method inference is unchanged to first order.
#' @param target Shrinkage target for \code{"eb"}: \code{"pooled"} (one
#'   accuracy for all classes, errors spread evenly), \code{"groups"}
#'   (separate rates for errors within and between \code{groups}, e.g. ISCO
#'   major groups), \code{"loglinear"} (Poisson quasi-independence model of
#'   the table with a diagonal and, with \code{groups}, a same-group term)
#'   or \code{"matrix"} (\code{target_matrix}). Ill-conditioned targets are
#'   refused.
#' @param groups Group of each category (in model order) for
#'   \code{target = "groups"} or \code{"loglinear"}.
#' @param target_matrix K x K column-stochastic target for
#'   \code{target = "matrix"}, e.g. a published confusion matrix.
#' @param concentration Prior concentration \eqn{c \ge 0} for
#'   \code{"eb"}, or \code{"ml"} (default) to maximise the
#'   Dirichlet-multinomial marginal likelihood.
#' @param alpha Pseudo-count per cell for \code{estimator = "dirichlet"}.
#' @param B Number of bootstrap replicates or posterior draws.
#' @param seed Optional seed for \code{variance = "bootstrap"} or
#'   \code{"posterior"}; the global random-number state is restored
#'   afterwards.
#' @param prevalence Source of the latent prevalence \eqn{\pi}:
#'   \code{"validation"} (proportions of the true category in the
#'   validation sample; Hajek or Horvitz--Thompson as set in
#'   \code{\link{validation_sample}}), \code{"em"} (maximum likelihood from
#'   the main study's proxy frequencies given \eqn{\hat\Pi}; only
#'   \eqn{\Pi} is transported from the validation sample) or
#'   \code{"inverse"} (\eqn{\hat\Pi^{-1}\hat p}; errors if outside the
#'   simplex).
#' @param variance \code{"delta"} (default) propagates the estimation of
#'   the misclassification probabilities into the regression covariance
#'   through the stacked sandwich; \code{"conditional"} treats them as
#'   known; \code{"bootstrap"} resamples the regression rows (keeping
#'   internal-audit membership) and an external audit (within strata) and
#'   refits everything \code{B} times; \code{"posterior"} (external audits,
#'   constant prevalence) draws \eqn{\Pi} from its Dirichlet posterior,
#'   refits with the draws treated as known and combines the conditional
#'   covariance and the between-draw spread as
#'   \eqn{\bar V + (1 + 1/B)\,\widehat{\mathrm{Var}}(\hat\psi^{(b)})}.
#' @param beta_equation Regression estimating equation for an internal
#'   validation sample. \code{"yi"}: validated rows contribute the score at
#'   their true category instead of the corrected score (Yi et al., 2019,
#'   eq. 23); valid when selection into the validation sample depends only
#'   on \eqn{(\hat Z, x)}. \code{"weighted"}: every row keeps the corrected
#'   score and validated rows add \eqn{d_i(S_i - U_i)}; valid for any known
#'   inclusion probabilities. \code{"auto"} (default) uses \code{"weighted"}
#'   when the validation sample has user-supplied design weights and
#'   \code{"yi"} otherwise.
#' @param kappa_max,sigma_min,min_class_n Diagnostic thresholds: the
#'   largest acceptable condition number and smallest acceptable singular
#'   value of \eqn{\hat\Pi}, and the smallest acceptable number of
#'   validation units per true category.
#' @param on_ill What to do when \code{\link{diagnose_mc}} flags a problem:
#'   \code{"warn"} (default), \code{"error"}, \code{"none"} or
#'   \code{"regularize"} (re-estimate \eqn{\Pi} with \code{estimator =
#'   "eb"}, target \code{"groups"} when \code{groups} is given and
#'   \code{"pooled"} otherwise, and warn).
#' @param em_maxit,em_tol Iteration limit and tolerance of the EM
#'   prevalence.
#' @param prevalence_model Optional one-sided formula for a covariate-
#'   dependent prevalence \eqn{\Pr(Z = \ell \mid x)}, a multinomial logit
#'   (binary logit for \eqn{K = 2}) fitted with \code{nnet::multinom}
#'   (Yi et al., 2019, Section 4.3). Use it when the latent category depends
#'   on covariates (e.g. occupation on region and industry): the constant
#'   prevalence of \code{"sub"}, \code{"ec"} and \code{"il"} then makes
#'   them inconsistent. The model is fitted on the validation sample
#'   (\code{prevalence = "validation"}; design-weighted) or by maximum
#'   likelihood from the main study's proxies given \eqn{\hat\Pi}
#'   (\code{prevalence = "em"}); \code{"il"} estimates it jointly.
#'   Variables are taken from the main-study data and, for an external
#'   audit with \code{prevalence = "validation"}, from
#'   \code{validation_sample(data = )}.
#' @return An object of class \code{"mc_control"}.
#' @seealso \code{\link{estimate_mc}}, \code{\link{mcglm}}
#' @export
control_mc <- function(estimator = c("mle", "eb", "dirichlet"),
                       prevalence = c("validation", "em", "inverse"),
                       variance = c("delta", "conditional", "bootstrap",
                                    "posterior"),
                       beta_equation = c("auto", "yi", "weighted"),
                       kappa_max = 100, sigma_min = 0.05, min_class_n = 10,
                       on_ill = c("warn", "error", "none", "regularize"),
                       em_maxit = 1000L, em_tol = 1e-10,
                       prevalence_model = NULL,
                       target = c("pooled", "groups", "loglinear", "matrix"),
                       groups = NULL, target_matrix = NULL,
                       concentration = "ml", alpha = 0.5,
                       B = 200L, seed = NULL) {
  estimator <- match.arg(estimator)
  target <- match.arg(target)
  if (!identical(concentration, "ml") &&
      (!is.numeric(concentration) || length(concentration) != 1L ||
       !is.finite(concentration) || concentration < 0))
    stop("concentration must be \"ml\" or one non-negative number.",
         call. = FALSE)
  if (target == "matrix" && estimator == "eb" && is.null(target_matrix))
    stop("target = 'matrix' needs target_matrix.", call. = FALSE)
  if (!is.numeric(B) || length(B) != 1L || B < 2)
    stop("B must be at least 2.", call. = FALSE)
  num_ok <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v) &&
    v >= 0
  if (!num_ok(kappa_max) || !num_ok(sigma_min) || !num_ok(min_class_n) ||
      !num_ok(em_maxit) || !num_ok(em_tol))
    stop("control_mc(): thresholds and EM settings must be single ",
         "non-negative numbers.", call. = FALSE)
  if (!is.null(prevalence_model) &&
      (!inherits(prevalence_model, "formula") || length(prevalence_model) != 2L))
    stop("prevalence_model must be a one-sided formula, e.g. ~ x1 + region.",
         call. = FALSE)
  prevalence <- match.arg(prevalence)
  if (!is.null(prevalence_model) && prevalence == "inverse")
    stop("prevalence = 'inverse' has no covariate-dependent version; use ",
         "'validation' or 'em' with prevalence_model.", call. = FALSE)
  structure(list(estimator = estimator,
                 prevalence = prevalence,
                 prevalence_model = prevalence_model,
                 variance = match.arg(variance),
                 beta_equation = match.arg(beta_equation),
                 kappa_max = kappa_max, sigma_min = sigma_min,
                 min_class_n = min_class_n, on_ill = match.arg(on_ill),
                 em_maxit = as.integer(em_maxit), em_tol = em_tol,
                 target = target, groups = groups,
                 target_matrix = target_matrix,
                 concentration = concentration, alpha = alpha,
                 B = as.integer(B), seed = seed),
            class = "mc_control")
}

#' Estimate misclassification probabilities from a validation sample
#'
#' Estimates the misclassification matrix
#' \eqn{\Pi_{j\ell} = \Pr(\hat Z = j \mid Z = \ell)}, the latent prevalence
#' \eqn{\pi} and the predictive matrix
#' \eqn{W_{j\ell} = \Pr(Z = \ell \mid \hat Z = j)} from a validation sample,
#' using its design weights, and returns their sandwich covariance. The
#' result can be passed to \code{\link{mcglm}} through
#' \code{mc(z, estimate)} or \code{validation = estimate}, in which case
#' the uncertainty is propagated into the regression estimates.
#'
#' @details
#' \eqn{\hat\Pi} uses design-weighted column ratios
#' \eqn{\sum_V d_i 1\{\hat Z_i = j, Z_i = \ell\} / \sum_V d_i 1\{Z_i = \ell\}},
#' which are the same under the Hajek and Horvitz--Thompson estimators.
#' With \code{prevalence = "validation"} the prevalence is
#' \eqn{\sum_V d_i 1\{Z_i = \ell\}} divided by \eqn{\sum_V d_i} (Hajek) or
#' by the frame size \eqn{N} (HT). For an audit stratified by the proxy
#' (\code{strata = "z_hat"}) with weights \eqn{N_j / n_j}, the implied
#' \eqn{W} is the within-stratum proportion of each true category.
#' Covariances come from the stacked estimating functions; a validation
#' sample with strata uses a stratum-centred, with-replacement meat for an
#' external sample, and the unit-level meat (exact for Poisson sampling,
#' conservative for fixed stratum sizes) for an internal one.
#'
#' @param validation A \code{\link{validation_sample}} or the list form
#'   \code{list(z, z_hat)} / \code{list(z, index)}.
#' @param z_hat Main-study proxies (factor, character or integer codes
#'   \code{0, ..., K-1}). Required for internal validation and for
#'   \code{prevalence = "em"} or \code{"inverse"}; a factor also fixes the
#'   category labels and their order.
#' @param K Number of categories (default: from the labels).
#' @param levels Optional category labels in model order. Defaults to
#'   \code{levels(z_hat)} for a factor (sorted values for character),
#'   \code{"0", ..., "K-1"} for integer codes, and otherwise the audit's
#'   own labels. Use it so that the estimate matches the levels of the
#'   \code{mc()} variable when \code{z_hat} is not given.
#' @param main_weights Optional main-study frequency weights.
#' @param control A \code{\link{control_mc}} object.
#' @param data Main-study data frame with the variables of
#'   \code{control$prevalence_model} (rows aligned with \code{z_hat}).
#' @return An object of class \code{"mc_estimate"} with components
#'   \code{Pi}, \code{pi}, \code{W} (labelled by category), \code{K},
#'   \code{levels}, \code{eta} (free parameters), \code{vcov} (their
#'   covariance), \code{diagnostics}, the bound validation sample and the
#'   control settings.
#' @references
#' Saerens, M., Latinne, P. and Decaestecker, C. (2002). Adjusting the
#' outputs of a classifier to new a priori probabilities: a simple
#' procedure. \emph{Neural Computation}, 14(1), 21--41.
#' \doi{10.1162/089976602753284446}.
#' @seealso \code{\link{validation_sample}}, \code{\link{diagnose_mc}},
#'   \code{\link{mcglm}}
#' @examples
#' set.seed(1)
#' z     <- rbinom(300, 1, 0.4)
#' z_hat <- ifelse(z == 1, rbinom(300, 1, 0.85), rbinom(300, 1, 0.1))
#' est <- estimate_mc(validation_sample(z, z_hat))
#' est
#' summary(est)
#' @export
estimate_mc <- function(validation, z_hat = NULL, K = NULL, levels = NULL,
                        main_weights = NULL, control = control_mc(),
                        data = NULL) {
  val <- as_validation_sample(validation)
  if (!inherits(control, "mc_control"))
    stop("control must be created by control_mc().", call. = FALSE)
  cats <- .mc_category_levels(val, z_hat, K, levels)
  .estimate_mc_codes(val, cats$z_hat, cats$levels, main_weights, control,
                     data = data)
}

#' estimate_mc() on 0-based main-study codes and known category labels
#' @keywords internal
.estimate_mc_codes <- function(val, z_hat, levels, main_weights = NULL,
                               control = control_mc(), data = NULL) {
  K <- length(levels)
  needs_main <- val$type == "internal" ||
    control$prevalence %in% c("em", "inverse")
  if (needs_main && is.null(z_hat))
    stop("estimate_mc() needs the main-study proxies z_hat for internal ",
         "validation and for prevalence = 'em' or 'inverse'.", call. = FALSE)
  if (!is.null(z_hat) && any(z_hat < 0L | z_hat >= K))
    stop("z_hat must be coded 0, ..., K-1 (K = ", K, ").", call. = FALSE)
  if (!is.null(main_weights) &&
      (length(main_weights) != length(z_hat) ||
       any(!is.finite(main_weights)) || any(main_weights <= 0)))
    stop("main_weights must be positive, one per main-study row.",
         call. = FALSE)

  vb   <- .mc_bind_validation(val, z_hat, K, main_weights, levels = levels)
  model <- control$prevalence_model
  w_main <- NULL
  w_v <- NULL
  if (is.null(model)) {
    nuis <- .mc_nuisance(vb, z_hat, K, control, main_weights)
    par <- nuis$map(nuis$eta)
  } else {
    des_m <- if (!is.null(data)) .mc_prev_design(model, data) else NULL
    if (!is.null(des_m) && !is.null(z_hat) && nrow(des_m$x) != length(z_hat))
      stop("estimate_mc(data = ) must have one row per main-study proxy.",
           call. = FALSE)
    if (is.null(des_m) && (needs_main || val$type == "internal"))
      .mc_prev_design(model, NULL)          # informative error
    w_main <- des_m$x
    w_v <- NULL
    des <- des_m
    if (control$prevalence == "validation") {
      if (val$type == "internal") {
        w_v <- w_main[vb$index, , drop = FALSE]
      } else {
        des_v <- .mc_prev_design(model, val$data, terms = des_m$terms,
                                 xlevels = des_m$xlevels)
        w_v <- des_v$x
        if (is.null(des)) des <- des_v
      }
    }
    nuis <- .mc_nuisance_model(vb, z_hat, K, control, main_weights, w_main,
                               w_v)
    nuis$prevalence$model <- list(formula = model, terms = des$terms,
                                  xlevels = des$xlevels,
                                  columns = colnames(des$x))
    par <- nuis$map(nuis$eta)
    probs <- .mc_prev_probs(par$alpha, if (is.null(w_main)) w_v else w_main)
    par$pi <- colMeans(probs)
    par$W <- NULL
  }
  dimnames(par$Pi) <- list(z_hat = levels, z = levels)
  if (!is.null(par$W)) dimnames(par$W) <- list(z_hat = levels, z = levels)
  names(par$pi)    <- levels
  est <- c(par,
           list(K = K, levels = levels, eta = nuis$eta, vcov = nuis$vcov,
                map = nuis$map, rows = nuis$rows, score = nuis$score,
                A = nuis$A, shrinkage = nuis$shrinkage,
                prevalence = nuis$prevalence, validation = vb,
                z_hat = if (needs_main) z_hat else NULL,
                w_main = w_main,
                w_v = if (is.null(model)) NULL else w_v,
                main_weights = main_weights, control = control,
                main_dependent = needs_main))
  class(est) <- "mc_estimate"
  est$diagnostics <- diagnose_mc(est)
  if (control$on_ill == "regularize" && length(est$diagnostics$problems) &&
      control$estimator == "mle") {
    reg <- control
    reg$estimator <- "eb"
    reg$target <- if (is.null(control$groups)) "pooled" else "groups"
    reg$on_ill <- "warn"
    warning("Estimated misclassification matrix: ",
            paste(est$diagnostics$problems, collapse = "; "),
            "; re-estimating with estimator = 'eb' (target = '", reg$target,
            "').", call. = FALSE)
    return(.estimate_mc_codes(val, z_hat, levels, main_weights, reg,
                              data = data))
  }
  .mc_on_ill(est$diagnostics, control$on_ill)
  est
}

#' Nuisance estimates, estimating functions and covariance
#' @keywords internal
.mc_nuisance <- function(vb, z_hat, K, control, wt = NULL) {
  s <- K - 1L
  n_pi <- K * s
  prevalence <- control$prevalence
  internal <- vb$type == "internal"
  n_m <- if (is.null(z_hat)) 0L else length(z_hat)
  wt_m <- if (is.null(wt)) rep(1, n_m) else wt
  d <- vb$d
  dw <- d * vb$w_v
  Zv <- outer(vb$z, 0:s, "==") * 1
  Hv <- outer(vb$proxy, 0:s, "==") * 1
  Hm <- if (n_m) outer(z_hat, 0:s, "==") * 1 else NULL
  N_tot <- vb$N
  if (is.null(N_tot) && !is.null(vb$N_h))
    N_tot <- sum(tapply(vb$N_h, vb$strata, `[`, 1L), na.rm = TRUE)

  # Pi uses design weights normalised to the number of validation units, so
  # that a shrinkage prior is weighed against the actual audit size.
  d_pi <- d * sum(vb$w_v) / sum(dw)
  C <- crossprod(Hv * (d_pi * vb$w_v), Zv)
  shrink <- .mc_shrink_setup(C, control)
  cc <- if (is.null(shrink)) 0 else shrink$c
  Tm <- if (is.null(shrink)) matrix(0, K, K) else shrink$T
  Pi_hat <- (C + cc * Tm) / rep(colSums(C) + cc, each = K)
  em_info <- NULL
  if (prevalence == "validation") {
    denom <- if (vb$estimator == "ht") N_tot else sum(dw)
    prev <- colSums(dw * Zv)[-1L] / denom
  } else if (prevalence == "em") {
    p_main <- colSums(wt_m * Hm) / sum(wt_m)
    em_info <- .mc_prevalence_em(p_main, Pi_hat, maxit = control$em_maxit,
                                 tol = control$em_tol)
    prev <- em_info$pi[-1L]
  } else {
    prev <- (colSums(wt_m * Hm) / sum(wt_m))[-1L]
  }
  eta <- c(as.numeric(Pi_hat[-1L, , drop = FALSE]), prev)
  q <- length(eta)

  map <- function(eta) {
    eta <- unname(eta)
    P  <- matrix(eta[seq_len(n_pi)], s, K)
    Pi <- rbind(1 - colSums(P), P)
    pr <- eta[n_pi + seq_len(s)]
    pi <- if (prevalence == "inverse") {
      tryCatch(as.numeric(solve(Pi, c(1 - sum(pr), pr))),
               error = function(e)
                 stop("prevalence = 'inverse': the estimated Pi is singular.",
                      call. = FALSE))
    } else {
      c(1 - sum(pr), pr)
    }
    list(Pi = Pi, pi = pi, W = .mcglm_predictive_matrix(Pi, pi))
  }

  # Per-unit estimating-function rows at eta: m (main units, or NULL when
  # the main study does not enter) and v (validation units).
  rows <- function(eta) {
    pr <- map(eta)
    v <- matrix(0, vb$n, q)
    for (l in seq_len(K))
      v[, (l - 1L) * s + seq_len(s)] <- d_pi * Zv[, l] *
        (Hv[, -1L, drop = FALSE] - rep(pr$Pi[-1L, l], each = vb$n))
    pc <- n_pi + seq_len(s)
    m <- NULL
    if (prevalence == "validation") {
      if (vb$estimator == "hajek") {
        v[, pc] <- d * (Zv[, -1L, drop = FALSE] - rep(pr$pi[-1L], each = vb$n))
      } else if (internal) {
        v[, pc] <- d * Zv[, -1L, drop = FALSE]
        m <- matrix(0, n_m, q)
        m[, pc] <- -rep(pr$pi[-1L], each = n_m)
      } else {
        v[, pc] <- d * Zv[, -1L, drop = FALSE] -
          rep(N_tot / vb$n * pr$pi[-1L], each = vb$n)
      }
    } else if (prevalence == "em") {
      m <- matrix(0, n_m, q)
      m[, pc] <- pr$W[z_hat + 1L, -1L, drop = FALSE] -
        rep(pr$pi[-1L], each = n_m)
    } else {
      m <- matrix(0, n_m, q)
      m[, pc] <- Hm[, -1L, drop = FALSE] - rep(unname(eta[pc]), each = n_m)
    }
    list(m = m, v = v)
  }
  total <- function(eta) {
    r <- rows(eta)
    out <- colSums(vb$w_v * r$v)
    if (!is.null(r$m)) out <- out + colSums(wt_m * r$m)
    out[seq_len(n_pi)] <- out[seq_len(n_pi)] +
      .mc_shrink_penalty(map(eta)$Pi, cc, Tm)
    out
  }
  A <- .mc_num_jacobian(total, eta)

  pi_hat <- map(eta)$pi
  if (prevalence == "inverse" &&
      (any(!is.finite(pi_hat)) || any(pi_hat <= 0 | pi_hat >= 1)))
    stop("prevalence = 'inverse': Pi^{-1} applied to the main-study proxy ",
         "frequencies gives prevalences (",
         paste(round(pi_hat, 4), collapse = ", "),
         ") outside (0, 1). Use prevalence = 'em'.", call. = FALSE)

  r <- rows(eta)
  B <- .mc_meat(r$m, r$v, vb, wt_m, n_m)
  V <- .mc_sandwich_or_na(A, B, "estimate_mc()")

  nms <- c(as.vector(outer(seq_len(s), 0:s,
                           function(j, l) sprintf("Pi[%d,%d]", j, l))),
           sprintf(if (prevalence == "inverse") "p_hat[%d]" else "pi[%d]",
                   seq_len(s)))
  names(eta) <- nms
  dimnames(V) <- list(nms, nms)
  list(eta = eta, map = map, rows = rows, score = total, A = A, vcov = V,
       shrinkage = shrink,
       prevalence = list(method = prevalence, em = em_info,
                         boundary = !is.null(em_info) && em_info$boundary))
}

#' Nuisance estimates with a covariate-dependent prevalence
#'
#' Parameters: the free entries of \eqn{\Pi} (as without a model) and the
#' multinomial-logit coefficients \eqn{\alpha} ((K - 1) x q, stored
#' column-major). Estimating functions for \eqn{\alpha}:
#' \itemize{
#'   \item \code{prevalence = "validation"}: validation units,
#'     \eqn{d_i\{1(Z_i = \ell) - \pi_\ell(x_i)\} x_i};
#'   \item \code{prevalence = "em"}: main units,
#'     \eqn{\{\Pr(Z = \ell \mid \hat Z_i, x_i) - \pi_\ell(x_i)\} x_i}, the
#'     score of \eqn{\log \sum_\ell \Pi_{\hat z_i \ell}\pi_\ell(x_i)}.
#' }
#' Starting values come from \code{nnet::multinom} (EM with multinom
#' M-steps for \code{"em"}); Newton steps then solve the stacked equations.
#' @keywords internal
.mc_nuisance_model <- function(vb, z_hat, K, control, wt, w_main, w_v) {
  s <- K - 1L
  n_pi <- K * s
  src <- control$prevalence
  n_m <- if (is.null(z_hat)) 0L else length(z_hat)
  wt_m <- if (is.null(wt)) rep(1, n_m) else wt
  d <- vb$d
  dw <- d * vb$w_v
  Zv <- outer(vb$z, 0:s, "==") * 1
  Hv <- outer(vb$proxy, 0:s, "==") * 1
  w_ref <- if (is.null(w_main)) w_v else w_main
  q_w <- ncol(w_ref)
  n_alpha <- s * q_w

  d_pi <- d * sum(vb$w_v) / sum(dw)
  C <- crossprod(Hv * (d_pi * vb$w_v), Zv)
  shrink <- .mc_shrink_setup(C, control)
  cc <- if (is.null(shrink)) 0 else shrink$c
  Tm <- if (is.null(shrink)) matrix(0, K, K) else shrink$T
  Pi_hat <- (C + cc * Tm) / rep(colSums(C) + cc, each = K)
  em_iter <- NULL
  if (src == "validation") {
    alpha <- .mc_multinom(vb$z, w_v, dw, K)
  } else {
    alpha <- .mc_multinom(z_hat, w_main, wt_m, K)
    rows_big <- rep(seq_len(n_m), K)
    ll_old <- -Inf
    for (em_iter in seq_len(min(control$em_maxit, 300L))) {
      joint <- Pi_hat[z_hat + 1L, , drop = FALSE] * .mc_prev_probs(alpha, w_main)
      ll <- sum(wt_m * log(rowSums(joint)))
      if (abs(ll - ll_old) < 1e-10 * (1 + abs(ll))) break
      ll_old <- ll
      post <- joint / rowSums(joint)
      alpha <- .mc_multinom(rep(0:s, each = n_m), w_main[rows_big, , drop = FALSE],
                            as.numeric(wt_m * post), K)
    }
  }

  map <- function(eta) {
    eta <- unname(eta)
    P  <- matrix(eta[seq_len(n_pi)], s, K)
    list(Pi = rbind(1 - colSums(P), P),
         alpha = matrix(eta[n_pi + seq_len(n_alpha)], s, q_w),
         pi = NULL, W = NULL)
  }
  rows <- function(eta) {
    pr <- map(eta)
    v <- matrix(0, vb$n, n_pi + n_alpha)
    for (l in seq_len(K))
      v[, (l - 1L) * s + seq_len(s)] <- d_pi * Zv[, l] *
        (Hv[, -1L, drop = FALSE] - rep(pr$Pi[-1L, l], each = vb$n))
    m <- NULL
    if (src == "validation") {
      R <- d * (Zv - .mc_prev_probs(pr$alpha, w_v))[, -1L, drop = FALSE]
      for (k in seq_len(q_w))
        v[, n_pi + (k - 1L) * s + seq_len(s)] <- R * w_v[, k]
    } else {
      P <- .mc_prev_probs(pr$alpha, w_main)
      joint <- pr$Pi[z_hat + 1L, , drop = FALSE] * P
      R <- (joint / rowSums(joint) - P)[, -1L, drop = FALSE]
      m <- matrix(0, n_m, n_pi + n_alpha)
      for (k in seq_len(q_w))
        m[, n_pi + (k - 1L) * s + seq_len(s)] <- R * w_main[, k]
    }
    list(m = m, v = v)
  }
  total <- function(eta) {
    r <- rows(eta)
    out <- colSums(vb$w_v * r$v)
    if (!is.null(r$m)) out <- out + colSums(wt_m * r$m)
    out[seq_len(n_pi)] <- out[seq_len(n_pi)] +
      .mc_shrink_penalty(map(eta)$Pi, cc, Tm)
    out
  }
  eta <- c(as.numeric(Pi_hat[-1L, , drop = FALSE]), as.numeric(alpha))
  scale <- sum(vb$w_v) + sum(wt_m)
  sol <- nleqslv::nleqslv(eta, function(e) total(e) / scale,
                          jac = function(e) .mc_num_jacobian(total, e) / scale,
                          control = list(maxit = 100, ftol = 1e-12))
  if (sol$termcd > 2)
    warning("estimate_mc(): the prevalence model did not converge ",
            "(termcd = ", sol$termcd, ").", call. = FALSE)
  eta <- sol$x
  A <- .mc_num_jacobian(total, eta)
  r <- rows(eta)
  B <- .mc_meat(r$m, r$v, vb, wt_m, n_m)
  V <- .mc_sandwich_or_na(A, B, "estimate_mc()")

  cols <- colnames(w_ref)
  if (is.null(cols)) cols <- paste0("w", seq_len(q_w))
  nms <- c(as.vector(outer(seq_len(s), 0:s,
                           function(j, l) sprintf("Pi[%d,%d]", j, l))),
           as.vector(outer(seq_len(s), cols,
                           function(l, k) sprintf("alpha[%d,%s]", l, k))))
  names(eta) <- nms
  dimnames(V) <- list(nms, nms)
  list(eta = eta, map = map, rows = rows, score = total, A = A, vcov = V,
       shrinkage = shrink,
       prevalence = list(method = src, em = list(iterations = em_iter),
                         boundary = FALSE))
}

#' Sandwich covariance, or NA with a warning when the bread is singular
#'
#' A singular bread arises at a boundary prevalence (EM) or with an empty
#' cell, where the delta-method covariance is not defined.
#' @keywords internal
.mc_sandwich_or_na <- function(A, B, label) {
  A_inv <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(A_inv)) {
    warning(label, ": the stacked estimating equations have a singular ",
            "Jacobian (e.g. a prevalence on the boundary); covariances are ",
            "set to NA.", call. = FALSE)
    return(matrix(NA_real_, nrow(B), ncol(B)))
  }
  A_inv %*% B %*% t(A_inv)
}

#' Central-difference Jacobian of a vector function
#' @keywords internal
.mc_num_jacobian <- function(f, x, h = 1e-6) {
  f0 <- f(x)
  J <- vapply(seq_along(x), function(k) {
    e <- replace(numeric(length(x)), k, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, numeric(length(f0)))
  matrix(J, length(f0), length(x))
}

#' Meat of a stacked sandwich over main and validation units
#'
#' @param G_m Rows of main-study units (n x k) or \code{NULL}.
#' @param G_v Rows of validation units (n_V x k).
#' @param vb Bound validation sample.
#' @param wt_m Main-study frequency weights.
#' @param n_m Number of main-study rows.
#' @details Internal validation rows are added to their main-study unit
#'   (one unit, one row) and the meat is the frequency-weighted sum of
#'   outer products. For an external sample the validation part is added
#'   separately: stratum-centred with factor \eqn{n_h/(n_h - 1)} when
#'   strata are given, the plain sum of outer products otherwise.
#' @keywords internal
.mc_meat <- function(G_m, G_v, vb, wt_m, n_m) {
  k <- ncol(G_v)
  if (vb$type == "internal") {
    G <- if (is.null(G_m)) matrix(0, n_m, k) else G_m
    G[vb$index, ] <- G[vb$index, , drop = FALSE] + G_v
    return(crossprod(G * wt_m, G))
  }
  B <- if (is.null(G_m)) matrix(0, k, k) else crossprod(G_m * wt_m, G_m)
  if (is.null(vb$strata)) return(B + crossprod(G_v * vb$w_v, G_v))
  for (h in levels(droplevels(vb$strata))) {
    Gh <- G_v[vb$strata == h, , drop = FALSE]
    nh <- nrow(Gh)
    if (nh > 1L) {
      Gc <- sweep(Gh, 2, colMeans(Gh))
      B <- B + nh / (nh - 1) * crossprod(Gc)
    } else {
      B <- B + crossprod(Gh)
    }
  }
  B
}

#' @export
print.mc_estimate <- function(x, digits = 4L, ...) {
  vb <- x$validation
  cat(sprintf("Estimated misclassification (K = %d) from an %s validation sample (n = %d)\n",
              x$K, vb$type, vb$n))
  cat(sprintf("Prevalence: %s; estimator: %s%s\n", x$prevalence$method,
              if (vb$estimator == "ht") "Horvitz-Thompson" else "Hajek",
              if (vb$user_weights) ", design weights" else ""))
  if (!is.null(x$shrinkage))
    cat(sprintf("Shrinkage: %s towards the %s target, concentration c = %.4g\n",
                x$shrinkage$estimator, x$shrinkage$target, x$shrinkage$c))
  cat("\nPi = P(z_hat = row | z = column):\n")
  print(round(x$Pi, digits))
  if (is.null(x$prevalence$model)) {
    cat("\nLatent prevalence pi:\n")
  } else {
    cat("\nPrevalence model ", deparse(x$prevalence$model$formula),
        "; average latent prevalence:\n", sep = "")
  }
  print(round(x$pi, digits))
  if (length(x$diagnostics$problems))
    cat("\nDiagnostics:", paste(x$diagnostics$problems, collapse = "; "), "\n")
  invisible(x)
}

#' @export
summary.mc_estimate <- function(object, ...) {
  se <- sqrt(pmax(diag(object$vcov), 0))
  tab <- cbind(Estimate = object$eta, `Std. Error` = se)
  structure(list(coefficients = tab, diagnostics = object$diagnostics,
                 estimate = object),
            class = "summary.mc_estimate")
}

#' @export
print.summary.mc_estimate <- function(x, digits = 4L, ...) {
  print(x$estimate, digits = digits)
  cat(sprintf("\nFree parameters (indices are 0-based codes of %s; row/category 0 implied):\n",
              paste(x$estimate$levels, collapse = ", ")))
  print(round(x$coefficients, digits))
  cat("\n")
  print(x$diagnostics)
  invisible(x)
}

#' @export
vcov.mc_estimate <- function(object, ...) object$vcov
