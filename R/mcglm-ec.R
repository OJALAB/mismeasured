# ---------------------------------------------------------------------------
# Expectation correction ("ec") and induced likelihood ("il")
#
# Reference: Yi, G. Y., Yan, Y., Liao, X. and Spiegelman, D. (2019).
# Parametric regression analysis with covariate misclassification in main
# study/validation study designs. Int. J. Biostat., 15(1), 20170002,
# Sections 3.1, 3.3 and 4.
#
# With class densities f(y | Z = l, x; psi) and prior weights
# P(Z = l | Z_hat_i), the expectation-corrected score is the posterior mean
# of the complete-data score,
#   U_i = sum_l omega_il S_il,  omega_il propto P(Z = l | Z_hat_i) f_il,
# which is also the score of the induced (mixture) log-likelihood
#   log sum_l P(Z = l | Z_hat_i) f_il                      (Appendix A).
# With known misclassification probabilities the two estimators coincide.
# They differ when the probabilities are estimated: "ec" plugs in estimates
# from the validation sample (two-stage), "il" maximises the joint
# likelihood of the main study and the validation sample over (psi, Pi, pi).
#
# The Jacobian of U follows Louis' identity:
#   dU_i/dtheta' = sum_l omega_il dS_il/dtheta' + sum_l omega_il S_il S_il'
#                  - U_i U_i'.
# Rows with a known category (internal validation) use omega = indicator,
# which gives the complete-data score and derivative.
# Gaussian models carry tau = log(sigma) as an extra parameter because the
# posterior weights depend on sigma.
# ---------------------------------------------------------------------------

#' Response, covariates and family for the ec / il class densities
#' @keywords internal
.mcglm_ec_model <- function(y, x, K, family) {
  family <- .normalize_family(family)
  if (!family$family %in% c("poisson", "binomial", "gaussian"))
    stop("ec/il support the poisson, binomial and gaussian families.",
         call. = FALSE)
  if (family$family == "poisson" && (any(y < 0) || any(y != floor(y))))
    stop("ec/il with family = poisson need a count response.", call. = FALSE)
  if (family$family == "binomial" && any(!y %in% c(0, 1)))
    stop("ec/il with family = binomial need a 0/1 response (grouped ",
         "binomial responses are not supported).", call. = FALSE)
  list(family = family, has_sigma = family$family == "gaussian", y = y,
       x = x, K = K, n = length(y), n_x = ncol(x), p = K - 1L + ncol(x),
       lfactorial_y = if (family$family == "poisson") lfactorial(y) else 0)
}

#' Restrict an ec / il model to a subset of rows
#' @keywords internal
.mcglm_ec_subset <- function(model, idx) {
  model$y <- model$y[idx]
  model$x <- model$x[idx, , drop = FALSE]
  model$n <- length(idx)
  if (length(model$lfactorial_y) > 1L) model$lfactorial_y <- model$lfactorial_y[idx]
  model
}

#' Per-class linear predictors, means and log densities
#' @keywords internal
.mcglm_ec_classes <- function(theta, model) {
  s <- model$K - 1L
  gamma <- c(0, theta[seq_len(s)])
  eta0 <- drop(model$x %*% theta[s + seq_len(model$n_x)])
  eta <- outer(eta0, gamma, "+")
  mu <- eta
  mu[] <- model$family$linkinv(eta)
  mud <- eta
  mud[] <- model$family$mu.eta(eta)
  s2 <- 1
  if (model$family$family == "poisson") {
    logf <- model$y * eta - mu - model$lfactorial_y
  } else if (model$family$family == "binomial") {
    log1pexp <- ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta)))
    logf <- model$y * eta - log1pexp
  } else {
    tau <- theta[model$p + 1L]
    s2 <- exp(2 * tau)
    logf <- -0.5 * log(2 * pi) - tau - (model$y - mu)^2 / (2 * s2)
  }
  list(eta = eta, mu = mu, mud = mud, logf = logf, s2 = s2)
}

#' Posterior-weighted score, log-likelihood and Louis Jacobian
#'
#' @param theta \code{c(psi, tau)} (tau only for gaussian).
#' @param model List from \code{.mcglm_ec_model}.
#' @param logprior n x K matrix of log prior class weights (ignored when
#'   \code{omega} is given).
#' @param omega Optional n x K matrix of fixed class weights (indicators for
#'   rows with a known category).
#' @param w Row weights for the Jacobian.
#' @param jac Compute the sum-scale Jacobian of \code{colSums(w * U)}.
#' @return List with \code{U} (n x q rows), \code{omega}, \code{loglik}
#'   (per row) and optionally \code{J}.
#' @keywords internal
.mcglm_ec_parts <- function(theta, model, logprior = NULL, omega = NULL,
                            w = NULL, jac = FALSE) {
  cl <- .mcglm_ec_classes(theta, model)
  K <- model$K
  s <- K - 1L
  n <- model$n
  p <- model$p
  q <- p + model$has_sigma
  if (is.null(omega)) {
    lj <- logprior + cl$logf
    m <- lj[cbind(seq_len(n), max.col(lj, ties.method = "first"))]
    ex <- exp(lj - m)
    tot <- rowSums(ex)
    omega <- ex / tot
    loglik <- m + log(tot)
  } else {
    loglik <- rowSums(ifelse(omega > 0, omega * cl$logf, 0))
  }
  phi <- cl$s2
  resid <- model$y - cl$mu
  xi_of <- function(l) {
    xi <- cbind(matrix(0, n, s), model$x)
    if (l > 1L) xi[, l - 1L] <- 1
    xi
  }
  U <- matrix(0, n, q)
  S_list <- vector("list", K)
  for (l in seq_len(K)) {
    S <- matrix(0, n, q)
    S[, seq_len(p)] <- xi_of(l) * (resid[, l] / phi)
    if (model$has_sigma) S[, q] <- resid[, l]^2 / phi - 1
    S_list[[l]] <- S
    U <- U + omega[, l] * S
  }
  out <- list(U = U, omega = omega, loglik = loglik)
  if (jac) {
    w <- if (is.null(w)) rep(1, n) else w
    J <- -crossprod(U * w, U)
    for (l in seq_len(K)) {
      S <- S_list[[l]]
      wl <- w * omega[, l]
      J <- J + crossprod(S * wl, S)
      xi <- xi_of(l)
      J[seq_len(p), seq_len(p)] <- J[seq_len(p), seq_len(p)] -
        crossprod(xi * (wl * cl$mud[, l] / phi), xi)
      if (model$has_sigma) {
        cross <- -2 * colSums(xi * (wl * resid[, l] / phi))
        J[seq_len(p), q] <- J[seq_len(p), q] + cross
        J[q, seq_len(p)] <- J[q, seq_len(p)] + cross
        J[q, q] <- J[q, q] - 2 * sum(wl * resid[, l]^2 / phi)
      }
    }
    out$J <- J
  }
  out
}

#' M-step for (psi, tau) given class weights
#'
#' Weighted GLM on the data expanded over the K classes.
#' @keywords internal
.mcglm_ec_mstep <- function(theta, model, omega, w) {
  K <- model$K
  s <- K - 1L
  xi_big <- do.call(rbind, lapply(seq_len(K), function(l) {
    xi <- cbind(matrix(0, model$n, s), model$x)
    if (l > 1L) xi[, l - 1L] <- 1
    xi
  }))
  wt_big <- as.numeric(w * omega)
  keep <- wt_big > 0
  fit <- suppressWarnings(stats::glm.fit(
    xi_big[keep, , drop = FALSE], rep(model$y, K)[keep], weights = wt_big[keep],
    start = theta[seq_len(model$p)], family = model$family,
    control = stats::glm.control(epsilon = 1e-10, maxit = 50)))
  psi <- fit$coefficients
  psi[is.na(psi)] <- theta[seq_len(model$p)][is.na(psi)]
  if (!model$has_sigma) return(psi)
  cl <- .mcglm_ec_classes(c(psi, theta[model$p + 1L]), model)
  s2 <- sum(w * omega * (model$y - cl$mu)^2) / sum(w)
  c(psi, 0.5 * log(s2))
}

#' Starting values (psi, tau) from a naive fit
#' @keywords internal
.mcglm_ec_start <- function(psi_naive, model, xi_hat, wt) {
  if (!model$has_sigma) return(psi_naive)
  w <- if (is.null(wt)) rep(1, model$n) else wt
  res <- model$y - drop(xi_hat %*% psi_naive)
  c(psi_naive, 0.5 * log(sum(w * res^2) / sum(w)))
}

#' Expectation correction with known prior class weights
#'
#' Solves \eqn{\sum_i w_i U_i(\theta) = 0} with \eqn{U_i} the posterior
#' mean of the complete-data score under prior weights
#' \eqn{\Pr(Z = \ell \mid \hat Z_i)}: EM warm start, then Newton with the
#' Louis Jacobian. Identical to the induced-likelihood estimator when the
#' misclassification probabilities are known.
#' @param W Predictive matrix \eqn{\Pr(Z = \ell \mid \hat Z = j)}.
#' @param Wi Optional n x K row-specific prior weights (overrides
#'   \code{W}; used under a prevalence model).
#' @return List with \code{coefficients} (psi), \code{theta}, \code{vcov}
#'   (psi block), \code{sigma}, convergence fields.
#' @keywords internal
.mcglm_fit_ec <- function(psi_naive, y, xi_hat, z_hat, x, K, family, W,
                          wt = NULL, label = "EC", Wi = NULL) {
  model <- .mcglm_ec_model(y, x, K, family)
  w <- if (is.null(wt)) rep(1, length(y)) else wt
  if (is.null(Wi)) Wi <- W[z_hat + 1L, , drop = FALSE]
  logprior <- log(Wi)
  theta <- .mcglm_ec_start(psi_naive, model, xi_hat, wt)
  sol <- .mcglm_ec_solve(theta, model, logprior, w)
  theta <- sol$x
  parts <- .mcglm_ec_parts(theta, model, logprior, w = w, jac = TRUE)
  V <- .mc_sandwich_or_na(parts$J, crossprod(parts$U * w, parts$U), label)
  if (sol$termcd > 2)
    warning(label, " solver did not converge (termcd = ", sol$termcd, ")")
  list(coefficients = theta[seq_len(model$p)], theta = theta,
       vcov = V[seq_len(model$p), seq_len(model$p), drop = FALSE],
       sigma = if (model$has_sigma) exp(theta[model$p + 1L]) else NULL,
       converged = sol$termcd <= 2, termcd = sol$termcd,
       iterations = sol$iter)
}

#' EM then Newton for the posterior-weighted score with fixed priors
#' @keywords internal
.mcglm_ec_solve <- function(theta, model, logprior, w, em_maxit = 50L) {
  N <- sum(w)
  f <- function(th) colSums(w * .mcglm_ec_parts(th, model, logprior)$U) / N
  jf <- function(th) .mcglm_ec_parts(th, model, logprior, w = w, jac = TRUE)$J / N
  newton <- function(th)
    nleqslv::nleqslv(th, f, jac = jf,
                     control = list(maxit = 200, ftol = 1e-11))
  em <- function(th, maxit) {
    ll_old <- -Inf
    for (it in seq_len(maxit)) {
      omega <- .mcglm_ec_parts(th, model, logprior)$omega
      th <- .mcglm_ec_mstep(th, model, omega, w)
      ll <- sum(w * .mcglm_ec_parts(th, model, logprior)$loglik)
      if (abs(ll - ll_old) < 1e-10 * (1 + abs(ll))) break
      ll_old <- ll
    }
    th
  }
  theta <- em(theta, em_maxit)
  sol <- newton(theta)
  if (sol$termcd > 2) sol <- newton(em(theta, 2000L))
  sol
}

#' Expectation correction with a validation sample (two-stage)
#'
#' Plugs the predictive matrix of an \code{\link{estimate_mc}} object into
#' the posterior-weighted score and solves it with
#' \code{.mcglm_fit_validated} (true-category score on internal validation
#' rows; stacked sandwich).
#' @keywords internal
.mcglm_fit_ec_validated <- function(psi_naive, y, xi_hat, z_hat, x, K,
                                    family, est, wt = NULL,
                                    control = control_mc()) {
  model <- .mcglm_ec_model(y, x, K, family)
  vb <- est$validation
  p <- model$p
  # Warm start: the known-probability solution at the estimated W.
  par0 <- .mc_with_rows(est$map(est$eta), z_hat, est$w_main)
  start <- .mcglm_fit_ec(psi_naive, y, xi_hat, z_hat, x, K, family, NULL,
                         wt = wt, Wi = par0$Wi)$theta
  U_fun <- function(theta, par)
    .mcglm_ec_parts(theta, model, log(par$Wi))$U
  J_fun <- function(theta, par, w)
    .mcglm_ec_parts(theta, model, log(par$Wi), w = w, jac = TRUE)$J
  true_fun <- NULL
  if (vb$type == "internal") {
    model_v <- .mcglm_ec_subset(model, vb$index)
    omega_v <- outer(vb$z, seq_len(K) - 1L, "==") * 1
    true_fun <- list(
      U = function(theta) .mcglm_ec_parts(theta, model_v, omega = omega_v)$U,
      J = function(theta, w)
        .mcglm_ec_parts(theta, model_v, omega = omega_v, w = w, jac = TRUE)$J)
  }
  fit <- .mcglm_fit_validated(start, est, y, x, z_hat, K, model$family, U_fun,
                              J_fun, wt = wt, control = control,
                              label = "EC", true_fun = true_fun)
  theta <- fit$coefficients
  fit$coefficients <- theta[seq_len(p)]
  fit$vcov <- fit$vcov[seq_len(p), seq_len(p), drop = FALSE]
  fit$sigma <- if (model$has_sigma) exp(theta[p + 1L]) else NULL
  fit$nuisance$sigma <- fit$sigma
  fit
}

#' Induced likelihood with a validation sample (joint estimation)
#'
#' Maximises the joint (pseudo-)likelihood of the main study and the
#' validation sample over \eqn{(\psi, \tau, \Pi, \pi)}:
#' \itemize{
#'   \item main-study rows: \eqn{\log \sum_\ell \Pi_{\hat z_i \ell}
#'     \pi_\ell f(y_i \mid \ell, x_i)};
#'   \item internal validation rows: the complete-data
#'     \eqn{\log f(y_i \mid z_i, x_i) + \log \Pi_{\hat z_i z_i} +
#'     \log \pi_{z_i}} (valid whenever selection into the audit depends
#'     only on observed \eqn{(y, \hat z, x)}; design weights are not
#'     needed);
#'   \item external validation units: \eqn{\tilde d_v \{\log
#'     \Pi_{\hat z_v z_v} + \log \pi_{z_v}\}}, without the prevalence term
#'     for \code{prevalence = "em"} (only \eqn{\Pi} transported), with
#'     design weights normalised to sum to \eqn{n_V}.
#' }
#' Fitted by EM (closed-form M-steps for \eqn{\Pi, \pi}) followed by Newton
#' steps; inference by the sandwich over units in the logit
#' parametrisation of \eqn{(\Pi, \pi)}.
#' @keywords internal
.mcglm_fit_il_validated <- function(psi_naive, y, xi_hat, z_hat, x, K,
                                    family, est, wt = NULL,
                                    control = control_mc()) {
  model  <- .mcglm_ec_model(y, x, K, family)
  .mcglm_check_estimate(est, z_hat, K, wt)
  vb <- est$validation
  n  <- model$n
  s  <- K - 1L
  p  <- model$p
  q  <- p + model$has_sigma
  w  <- if (is.null(wt)) rep(1, n) else wt
  internal <- vb$type == "internal"
  prev_mode <- est$prevalence$method
  if (!internal && prev_mode == "inverse")
    stop("method = 'il' estimates the prevalence by maximum likelihood; ",
         "use control_mc(prevalence = 'validation') or 'em'.", call. = FALSE)
  v_prev <- !internal && prev_mode == "validation"
  has_model <- !is.null(est$prevalence$model)
  w_main <- est$w_main
  w_v <- if (has_model && v_prev) est$w_v else NULL
  q_w <- if (has_model) ncol(w_main) else 1L

  Hm <- outer(z_hat, 0:s, "==") * 1
  omega_fix <- NULL
  if (internal) {
    omega_fix <- matrix(NA_real_, n, K)
    omega_fix[vb$index, ] <- outer(vb$z, 0:s, "==") * 1
  }
  if (!internal) {
    dt <- vb$d * vb$n / sum(vb$d)
    Zv <- outer(vb$z, 0:s, "==") * 1
    Hv <- outer(vb$proxy, 0:s, "==") * 1
  }

  n_a <- s * K                       # Pi logits
  n_b <- s * q_w                     # prevalence logits or model coefficients
  b_cols <- q + n_a + seq_len(n_b)
  softmax_rows <- function(eta) {
    eta <- eta - apply(eta, 1L, max)
    ex <- exp(eta)
    ex / rowSums(ex)
  }
  unpack <- function(th) {
    a <- matrix(th[q + seq_len(n_a)], s, K)
    Pi <- apply(rbind(0, a), 2, function(v) exp(v - max(v)) / sum(exp(v - max(v))))
    b <- th[b_cols]
    if (has_model) {
      alpha <- matrix(b, s, q_w)
      list(theta = th[seq_len(q)], Pi = Pi, alpha = alpha,
           P = .mc_prev_probs(alpha, w_main),
           P_v = if (v_prev) .mc_prev_probs(alpha, w_v) else NULL)
    } else {
      pi <- softmax_rows(matrix(c(0, b), 1L))[1L, ]
      list(theta = th[seq_len(q)], Pi = Pi, pi = pi,
           P = matrix(pi, n, K, byrow = TRUE),
           P_v = if (v_prev) matrix(pi, vb$n, K, byrow = TRUE) else NULL)
    }
  }
  Pi_logits <- function(Pi) {
    Pi <- pmax(Pi, 1e-8)
    Pi <- sweep(Pi, 2, colSums(Pi), "/")
    as.numeric(log(Pi[-1L, , drop = FALSE]) - rep(log(Pi[1L, ]), each = s))
  }
  pi_logits <- function(pi) {
    pi <- pmax(pi, 1e-8)
    log(pi[-1L] / sum(pi)) - log(pi[1L] / sum(pi))
  }
  # Score rows of the prevalence parameters from residuals R (rows x s).
  prev_scores <- function(R, wmat) {
    if (!has_model) return(R)
    out <- matrix(0, nrow(R), n_b)
    for (k in seq_len(q_w)) out[, (k - 1L) * s + seq_len(s)] <- R * wmat[, k]
    out
  }

  # Main-study rows (internal validation rows use their true category).
  main_parts <- function(th, jac = FALSE) {
    u <- unpack(th)
    logprior <- log(u$Pi)[z_hat + 1L, , drop = FALSE] + log(u$P)
    pp <- .mcglm_ec_parts(u$theta, model, logprior, w = w, jac = jac)
    if (internal) {
      model_v <- .mcglm_ec_subset(model, vb$index)
      pv <- .mcglm_ec_parts(u$theta, model_v, omega = omega_fix[vb$index, ],
                            w = w[vb$index], jac = jac)
      pp$U[vb$index, ] <- pv$U
      pp$omega[vb$index, ] <- omega_fix[vb$index, ]
      pp$loglik[vb$index] <- pv$loglik +
        log(u$Pi[cbind(z_hat[vb$index] + 1L, vb$z + 1L)]) +
        log(u$P[cbind(vb$index, vb$z + 1L)])
      if (jac) {
        pn <- .mcglm_ec_parts(u$theta, model_v, logprior[vb$index, , drop = FALSE],
                              w = w[vb$index], jac = TRUE)
        pp$J <- pp$J - pn$J + pv$J
      }
    }
    G <- matrix(0, n, q + n_a + n_b)
    G[, seq_len(q)] <- pp$U
    om <- pp$omega
    for (l in seq_len(K))
      G[, q + (l - 1L) * s + seq_len(s)] <- om[, l] *
        (Hm[, -1L, drop = FALSE] - rep(u$Pi[-1L, l], each = n))
    G[, b_cols] <- prev_scores((om - u$P)[, -1L, drop = FALSE], w_main)
    list(G = G, loglik = pp$loglik, omega = om, J = pp$J, u = u)
  }
  ext_rows <- function(th) {
    if (internal) return(NULL)
    u <- unpack(th)
    G <- matrix(0, vb$n, q + n_a + n_b)
    for (l in seq_len(K))
      G[, q + (l - 1L) * s + seq_len(s)] <- dt * Zv[, l] *
        (Hv[, -1L, drop = FALSE] - rep(u$Pi[-1L, l], each = vb$n))
    if (v_prev)
      G[, b_cols] <- prev_scores(dt * (Zv - u$P_v)[, -1L, drop = FALSE], w_v)
    G
  }
  total_score <- function(th) {
    out <- colSums(w * main_parts(th)$G)
    if (!internal) out <- out + colSums(ext_rows(th))
    out
  }
  total_loglik <- function(mp, u) {
    ll <- sum(w * mp$loglik)
    if (!internal) {
      ll <- ll + sum(dt * log(u$Pi[cbind(vb$proxy + 1L, vb$z + 1L)]))
      if (v_prev) ll <- ll + sum(dt * log(u$P_v[cbind(seq_len(vb$n), vb$z + 1L)]))
    }
    ll
  }

  # EM: closed-form M-step for Pi, closed form (constant) or weighted
  # nnet::multinom (prevalence model) for the prevalence, weighted GLM for
  # (psi, tau).
  par_est <- est$map(est$eta)
  th <- c(.mcglm_ec_start(psi_naive, model, xi_hat, wt), Pi_logits(unname(est$Pi)),
          if (has_model) as.numeric(par_est$alpha) else pi_logits(unname(est$pi)))
  ll_old <- -Inf
  for (it in seq_len(500L)) {
    mp <- main_parts(th)
    om <- mp$omega
    cnt <- crossprod(Hm * w, om)            # K x K: [j, l]
    if (!internal) cnt <- cnt + crossprod(Hv * dt, Zv)
    Pi_new <- sweep(cnt, 2, colSums(cnt), "/")
    if (has_model) {
      codes <- rep(0:s, each = n)
      design <- w_main[rep(seq_len(n), K), , drop = FALSE]
      wts <- as.numeric(w * om)
      if (v_prev) {
        codes <- c(codes, vb$z)
        design <- rbind(design, w_v)
        wts <- c(wts, dt)
      }
      b_new <- as.numeric(.mc_multinom(codes, design, wts, K))
    } else {
      pc <- colSums(w * om)
      if (v_prev) pc <- pc + colSums(dt * Zv)
      b_new <- pi_logits(pc / sum(pc))
    }
    th_psi <- .mcglm_ec_mstep(th[seq_len(q)], model, om, w)
    th <- c(th_psi, Pi_logits(Pi_new), b_new)
    ll <- total_loglik(main_parts(th), unpack(th))
    if (abs(ll - ll_old) < 1e-11 * (1 + abs(ll))) break
    ll_old <- ll
  }
  sol <- nleqslv::nleqslv(th, function(t) total_score(t) / sum(w),
                          jac = function(t)
                            .mc_num_jacobian(total_score, t) / sum(w),
                          control = list(maxit = 100, ftol = 1e-10))
  if (sol$termcd > 2)
    warning("IL solver did not converge (termcd = ", sol$termcd, ")")
  th <- sol$x

  A <- .mc_num_jacobian(total_score, th)
  G_m <- main_parts(th)$G
  G_v <- if (internal) matrix(0, vb$n, ncol(G_m)) else ext_rows(th)
  vb_meat <- vb
  if (!internal) vb_meat$w_v <- rep(1, vb$n)
  B <- .mc_meat(G_m, G_v, vb_meat, w, n)
  Q <- length(th)
  if (control$variance == "conditional") {
    V <- matrix(NA_real_, Q, Q)
    V[seq_len(q), seq_len(q)] <- .mc_sandwich_or_na(
      A[seq_len(q), seq_len(q), drop = FALSE],
      B[seq_len(q), seq_len(q), drop = FALSE], "IL")
  } else {
    V <- .mc_sandwich_or_na(A, B, "IL")
  }

  u <- unpack(th)
  # Delta method from the logits to the reported nuisance parameters (free
  # entries of Pi, then pi or the prevalence-model coefficients).
  to_eta <- function(lg) {
    uu <- unpack(c(th[seq_len(q)], lg))
    c(as.numeric(uu$Pi[-1L, , drop = FALSE]),
      if (has_model) as.numeric(uu$alpha) else uu$pi[-1L])
  }
  Jp <- .mc_num_jacobian(to_eta, th[-seq_len(q)])
  V_eta <- Jp %*% V[-seq_len(q), -seq_len(q), drop = FALSE] %*% t(Jp)
  eta_names <- names(est$eta)
  if (!has_model && est$prevalence$method == "inverse")
    eta_names <- NULL
  if (is.null(eta_names) || length(eta_names) != nrow(V_eta))
    eta_names <- c(as.vector(outer(seq_len(s), 0:s,
                                   function(j, l) sprintf("Pi[%d,%d]", j, l))),
                   sprintf("pi[%d]", seq_len(s)))
  dimnames(V_eta) <- list(eta_names, eta_names)
  lev <- est$levels
  Pi <- u$Pi
  dimnames(Pi) <- list(z_hat = lev, z = lev)
  pi <- stats::setNames(if (has_model) colMeans(u$P) else u$pi, lev)
  W <- NULL
  if (!has_model) {
    W <- .mcglm_predictive_matrix(u$Pi, u$pi)
    dimnames(W) <- list(z_hat = lev, z = lev)
  }
  list(coefficients = th[seq_len(p)],
       vcov = V[seq_len(p), seq_len(p), drop = FALSE],
       sigma = if (model$has_sigma) exp(th[p + 1L]) else NULL,
       converged = sol$termcd <= 2, termcd = sol$termcd,
       iterations = sol$iter, em_iterations = it,
       nuisance = list(Pi = Pi, pi_z = pi, W = W,
                       alpha = if (has_model) u$alpha else NULL,
                       vcov = V_eta,
                       sigma = if (model$has_sigma) exp(th[p + 1L]) else NULL,
                       prevalence = if (internal) "joint" else prev_mode,
                       variance = control$variance))
}
