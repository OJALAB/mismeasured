# ---------------------------------------------------------------------------
# Replication-based variances for estimated misclassification probabilities
# (control_mc(variance = "bootstrap" / "posterior")).
#
# bootstrap: resample the regression rows with replacement (an internal
#   audit's membership travels with its rows, i.e. Poisson sampling of the
#   audit within the main study) and an external audit within its strata;
#   re-estimate the probabilities and refit every corrected method.
# posterior: for an external audit with a constant prevalence, draw each
#   column of Pi from its Dirichlet posterior (prior c T under shrinkage,
#   Jeffreys 0.5 otherwise) and the prevalence likewise (or by EM / inversion
#   from the main study), refit with the draws treated as known, and combine
#   V = mean_b V_b + (1 + 1/B) Var_b(psi_b)   (Rubin's rule).
# ---------------------------------------------------------------------------

#' Evaluate an expression with a temporary seed
#' @keywords internal
.mc_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  force(expr)
}

#' Resample a validation sample along with the main study
#'
#' @param val An \code{"mc_validation"} object.
#' @param idx Resampled main-study row numbers.
#' @return The validation sample of the resampled data.
#' @keywords internal
.mc_resample_validation <- function(val, idx) {
  if (val$type == "internal") {
    pos <- which(idx %in% val$index)
    src <- match(idx[pos], val$index)
    return(validation_sample(val$z[src], index = pos,
                             weights = val$weights[src],
                             strata = if (val$strata_type == "z_hat") "z_hat"
                                      else val$strata[src],
                             estimator = val$estimator, N = NULL,
                             N_strata = val$N_strata))
  }
  groups <- switch(val$strata_type, z_hat = val$z_hat, user = val$strata,
                   none = rep(1L, val$n))
  pick <- unlist(lapply(split(seq_len(val$n), groups), function(g)
    g[sample.int(length(g), length(g), replace = TRUE)]), use.names = FALSE)
  validation_sample(val$z[pick], val$z_hat[pick], weights = val$weights[pick],
                    strata = if (val$strata_type == "z_hat") "z_hat"
                             else val$strata[pick],
                    estimator = val$estimator, N = val$N,
                    N_strata = val$N_strata,
                    data = if (is.null(val$data)) NULL else
                      val$data[pick, , drop = FALSE])
}

#' Dirichlet draw
#' @keywords internal
.mc_rdirichlet <- function(a) {
  g <- stats::rgamma(length(a), shape = a)
  g / sum(g)
}

#' Replace covariances by bootstrap or posterior replication
#'
#' @param out The fitted \code{"mcglm"} object.
#' @param args The evaluated arguments of the original \code{mcglm()} call.
#' @return \code{out} with \code{vcov}, \code{se} and \code{replicates}
#'   updated for the corrected methods.
#' @keywords internal
.mcglm_replicate <- function(out, args) {
  ctrl <- args$mc_control
  type <- ctrl$variance
  B <- ctrl$B
  est <- out$mc_estimate
  methods <- setdiff(names(out$coefficients), "naive")
  supplied <- !vapply(args[c("Pi", "pi_z", "p01", "p10", "c1", "c2")],
                      is.null, logical(1))
  if (!is.null(out$formula) &&
      !is.null(.mcglm_parse_formula(out$formula, args$data,
                                    environment(out$formula))$Pi))
    supplied["Pi"] <- TRUE
  if (any(supplied))
    stop("variance = '", type, "' re-estimates the probabilities; do not ",
         "supply ", paste(names(supplied)[supplied], collapse = ", "),
         " (use mc(z) without a matrix).", call. = FALSE)
  if (type == "posterior") {
    if (est$validation$type != "external")
      stop("variance = 'posterior' is for external audits; use 'bootstrap' ",
           "or 'delta' for an internal one.", call. = FALSE)
    if (!is.null(est$prevalence$model))
      stop("variance = 'posterior' needs a constant prevalence; use ",
           "'bootstrap' or 'delta' with a prevalence model.", call. = FALSE)
  }
  if (inherits(args$validation, "mc_estimate") || is.null(args$validation))
    stop("variance = '", type, "' needs the validation sample itself ",
         "(validation_sample() or a list), not an estimate_mc() object or ",
         "mc(z, estimate).", call. = FALSE)
  val <- as_validation_sample(args$validation)
  formula_fit <- !is.null(out$formula)
  n <- out$n
  base <- args
  base$mc_control$variance <- "conditional"
  base$mc_control$on_ill <- "none"
  refit <- function(a) tryCatch(suppressWarnings(do.call(mcglm, a)),
                                error = function(e) NULL)
  keep_methods <- methods

  draws <- .mc_with_seed(ctrl$seed, lapply(seq_len(B), function(b) {
    a <- base
    if (type == "bootstrap") {
      idx <- sample.int(n, n, replace = TRUE)
      if (formula_fit) {
        a$data <- args$data[idx, , drop = FALSE]
      } else {
        a$formula <- args$formula[idx]
        a$z_hat <- args$z_hat[idx]
        a$x <- as.matrix(args$x)[idx, , drop = FALSE]
      }
      if (!is.null(args$weights)) a$weights <- args$weights[idx]
      a$validation <- .mc_resample_validation(val, idx)
    } else {
      pr <- .mc_posterior_draw(est, out)
      a$validation <- NULL
      a$Pi <- pr$Pi
      a$pi_z <- if (out$K == 2L) pr$pi[2L] else pr$pi
      a$method <- setdiff(args$method, "il")
    }
    fit <- refit(a)
    if (is.null(fit)) return(NULL)
    list(coef = fit$coefficients, vcov = fit$vcov)
  }))
  ok <- !vapply(draws, is.null, logical(1))
  if (sum(ok) < 2L)
    stop("variance = '", type, "': fewer than two replicates could be ",
         "fitted.", call. = FALSE)
  draws <- draws[ok]
  # il estimates Pi jointly (a posterior draw has no meaning for it); the
  # conditional covariances of bca/bcm are drifting-regime formulas that
  # understate the fixed-misclassification variance, so Rubin's rule would
  # too. They keep their delta-method covariance.
  if (type == "posterior")
    keep_methods <- setdiff(methods, c("il", "bca", "bcm"))
  reps <- list()
  for (m in keep_methods) {
    coefs <- do.call(rbind, lapply(draws, function(d) d$coef[[m]]))
    reps[[m]] <- coefs
    V <- stats::cov(coefs)
    if (type == "posterior") {
      Vb <- Reduce(`+`, lapply(draws, function(d) d$vcov[[m]])) / length(draws)
      V <- Vb + (1 + 1 / length(draws)) * V
    }
    dimnames(V) <- dimnames(out$vcov[[m]])
    out$vcov[[m]] <- V
    out$se[[m]] <- stats::setNames(sqrt(pmax(diag(V), 0)), rownames(V))
  }
  out$replicates <- list(type = type, B = B, failed = sum(!ok),
                         coefficients = reps,
                         delta_kept = setdiff(methods, keep_methods))
  out
}

#' One draw of (Pi, pi) from the posterior of an external audit
#' @keywords internal
.mc_posterior_draw <- function(est, out) {
  vb <- est$validation
  K <- est$K
  s <- K - 1L
  Zv <- outer(vb$z, 0:s, "==") * 1
  Hv <- outer(vb$proxy, 0:s, "==") * 1
  dw <- vb$d * vb$w_v
  d_pi <- vb$d * sum(vb$w_v) / sum(dw)
  C <- crossprod(Hv * (d_pi * vb$w_v), Zv)
  prior <- if (is.null(est$shrinkage)) matrix(0.5, K, K) else
    est$shrinkage$c * est$shrinkage$T
  Pi <- vapply(seq_len(K), function(l) .mc_rdirichlet(C[, l] + prior[, l]),
               numeric(K))
  method <- est$prevalence$method
  if (method == "validation") {
    counts <- colSums(Zv * (vb$d * vb$w_v)) * vb$n / sum(dw)
    pi <- .mc_rdirichlet(counts + 0.5)
  } else {
    zh <- out$z_hat
    w <- if (is.null(out$weights)) rep(1, length(zh)) else out$weights
    p <- as.numeric(tapply(w, factor(zh, levels = 0:s), sum))
    p[is.na(p)] <- 0
    p <- p / sum(p)
    pi <- if (method == "em") .mc_prevalence_em(p, Pi)$pi else
      as.numeric(solve(Pi, p))
  }
  list(Pi = Pi, pi = pi)
}
