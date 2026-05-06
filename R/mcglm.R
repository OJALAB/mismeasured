#' Bias-Corrected GLM Estimation with Misclassified Covariates
#'
#' Fits a GLM in which one covariate is observed through a misclassified
#' proxy and applies the bias-correction estimators of Battaglia,
#' Christensen, Hansen and Sacher (2025). Five estimators are available:
#' \emph{naive} (uncorrected proxy-score), \emph{BCA} (additive bias
#' correction), \emph{BCM} (multiplicative bias correction), \emph{CS}
#' (corrected score), and \emph{one-step} (joint mixture likelihood via
#' \pkg{RTMB}). Both binary and multicategory latent regressors are
#' supported.
#'
#' @section Model and notation:
#' Let \eqn{Y_i} denote a scalar response, \eqn{x_i \in \mathbb{R}^r} a
#' vector of cleanly observed covariates (typically including an
#' intercept), and \eqn{Z_i \in \{0, 1, \dots, K-1\}} a latent categorical
#' regressor that is observed only through a proxy \eqn{\hat Z_i}.
#' For the binary case (\eqn{K = 2}) we let
#' \eqn{\xi_i = (Z_i, x_i^\top)^\top} and
#' \eqn{\hat\xi_i = (\hat Z_i, x_i^\top)^\top}; for \eqn{K > 2},
#' \eqn{Z_i} is dummy-encoded with baseline category 0 and
#' \eqn{\xi_i = (d_i^\top, x_i^\top)^\top}.
#' The parameter is \eqn{\psi = (\gamma^\top, \alpha^\top)^\top}, where
#' \eqn{\gamma \in \mathbb{R}^{K-1}} are the coefficients on the
#' (dummy-coded) latent regressor and \eqn{\alpha \in \mathbb{R}^r}
#' the coefficients on \eqn{x_i}. The conditional mean is
#' \eqn{E[Y_i | Z_i, x_i] = \mu(\psi^\top \xi_i)} for
#' \eqn{\mu = g^{-1}} the inverse canonical link.
#'
#' \emph{Misclassification mechanism.}
#' For \eqn{K = 2}, the mechanism is summarized by the false-positive
#' rate \eqn{p_{01} = \Pr(\hat Z = 1 \mid Z = 0)}, the false-negative
#' rate \eqn{p_{10} = \Pr(\hat Z = 0 \mid Z = 1)} and the prevalence
#' \eqn{\pi = \Pr(Z = 1)}. For \eqn{K > 2}, by the \eqn{K \times K}
#' matrix \eqn{\Pi_{j\ell} = \Pr(\hat Z = j - 1 \mid Z = \ell - 1)}
#' (column-stochastic) and the prevalence vector
#' \eqn{(\pi_0, \dots, \pi_{K-1})}. Misclassification is assumed
#' nondifferential: \eqn{\hat Z \perp\!\!\!\perp (Y, x) \mid Z}.
#'
#' \emph{Estimators.}
#' Let \eqn{\hat\psi} solve the proxy score
#' \eqn{\hat U_n(\psi) = n^{-1} \sum_i \hat\xi_i \{Y_i - \mu(\psi^\top \hat\xi_i)\} = 0}
#' (this is the naive GLM estimator). Define the per-observation drift
#' \eqn{m_i(\psi) = (-c_1 \delta_i(\psi),\, -c_2 \delta_i(\psi)\, x_i)^\top}
#' with \eqn{\delta_i(\psi) = \mu(\gamma + \alpha^\top x_i) - \mu(\alpha^\top x_i)},
#' \eqn{c_1 = p_{01}(1 - \pi)} and
#' \eqn{c_2 = p_{01}(1 - \pi) - p_{10} \pi}, and let
#' \eqn{\hat I(\psi) = n^{-1} \sum_i \dot\mu(\psi^\top \hat\xi_i) \hat\xi_i \hat\xi_i^\top}
#' and \eqn{\hat M(\psi) = \partial \hat m / \partial \psi^\top}. Then:
#'
#' \describe{
#'   \item{BCA:}{\eqn{\hat\psi_{\mathrm{bca}} = \hat\psi - \hat I(\hat\psi)^{-1} \hat m(\hat\psi)}}
#'   \item{BCM:}{\eqn{\hat\psi_{\mathrm{bcm}} = \hat\psi - (\hat I(\hat\psi) + \hat M(\hat\psi))^{-1} \hat m(\hat\psi)}}
#'   \item{CS:}{solves \eqn{n^{-1} \sum_i \phi_i(\psi) = 0} with
#'     \eqn{\phi_i(\psi) = \hat\xi_i \{Y_i - \mu(\psi^\top \hat\xi_i)\} - m_i(\psi)}}
#'   \item{one-step:}{maximizes the integrated mixture likelihood
#'     \eqn{\prod_i \sum_\ell \pi_\ell\, \Pi_{\hat z_i+1, \ell+1}\, f(Y_i \mid Z = \ell, \psi)}
#'     by automatic differentiation (\pkg{RTMB}).}
#' }
#'
#' Multicategory analogues replace \eqn{(p_{01}, p_{10}, \pi)} with
#' \eqn{(\Pi, \pi_z)} throughout; see Appendix B of Battaglia et al.
#' (2025) for explicit formulas.
#'
#' \emph{Inference.}
#' All five estimators report an asymptotic variance through
#' \code{\link{vcov.mcglm}}: \eqn{A^{-1} C A^{-1}} for naive / BCA / BCM
#' (Theorems on the two-step expansion and on bias-corrected estimators
#' under drifting misclassification), \eqn{J^{-1} S J^{-\top}} for CS
#' (Z-estimator sandwich, with \eqn{J = -(\hat I + \hat M)},
#' \eqn{S = E[\phi_i \phi_i^\top]}), and the inverse Hessian of the
#' integrated log-likelihood for one-step. See \code{vcov_corrected}
#' below for an alternative, more conservative BCA/BCM sandwich.
#'
#' @section Required inputs by method:
#' Only the chosen \code{method}s' inputs are validated; you do not have
#' to supply parameters for methods you are not fitting.
#'
#' \tabular{lll}{
#'   \strong{Method}    \tab \strong{Binary (K = 2)}                    \tab \strong{Multicategory (K > 2)} \cr
#'   \code{"naive"}    \tab \emph{none}                                 \tab \emph{none} \cr
#'   \code{"bca"}, \code{"bcm"}, \code{"cs"} \tab \code{c1} and \code{c2}, or \code{p01}/\code{p10}/\code{pi_z}, or \code{Pi} \tab \code{Pi} and \code{pi_z} \cr
#'   \code{"onestep"} (\code{fix_omega = FALSE}) \tab \emph{none} (mixture weights estimated) \tab \emph{none} \cr
#'   \code{"onestep"} (\code{fix_omega = TRUE})  \tab \code{p01}, \code{p10}, \code{pi_z} \tab \code{Pi}, \code{pi_z}
#' }
#'
#' \emph{Auto-derivations performed by} \code{mcglm}\emph{:}
#' \itemize{
#'   \item If \code{Pi} is supplied with \code{nrow(Pi) == 2},
#'     \code{p01 = Pi[2, 1]} and \code{p10 = Pi[1, 2]} are extracted.
#'   \item If \code{Pi} and \code{z_hat} are available but \code{pi_z}
#'     is not, \code{pi_z} is estimated by Bayesian inversion of the
#'     observed proxy frequencies (\code{.mcglm_estimate_pi_z}).
#'   \item If \code{p01}, \code{p10}, \code{pi_z} are available but
#'     \code{c1}, \code{c2} are not, the latter are computed from the
#'     formulas above.
#' }
#' Supplying \code{c1} and \code{c2} directly is an "advanced" path
#' that bypasses these checks. The formula interface
#' (\code{y ~ mc(z, Pi) + ...}) extracts \code{Pi} automatically and is
#' the recommended entry point.
#'
#' @param formula Either a formula of the form
#'   \code{y ~ mc(z, Pi) + x1 + x2}, where \code{\link{mc}()} marks the
#'   misclassified covariate \eqn{\hat Z_i} together with its
#'   misclassification matrix \eqn{\Pi}; or, when using the matrix
#'   interface, a numeric response vector of length \eqn{n}.
#' @param data A data frame (required for the formula interface). The
#'   variable named in \code{mc()} should be a \code{\link{factor}} or
#'   integer-valued vector with levels coded as \eqn{0, 1, \dots, K-1}.
#' @param family A \code{\link[stats]{family}} object or one of
#'   the strings \code{"poisson"}, \code{"binomial"}, \code{"gaussian"}
#'   (any \eqn{K}, all five methods), or \code{"multinomial"}
#'   (\code{naive} and \code{onestep} only).
#' @param method Character vector of estimators to fit. Any subset of
#'   \code{c("naive", "bca", "bcm", "cs", "onestep")}; the default is
#'   the four analytical estimators.
#' @param p01 False-positive rate \eqn{p_{01} = \Pr(\hat Z = 1 \mid Z = 0)}
#'   (\code{K = 2} only). Auto-extracted from \code{Pi} when supplied.
#' @param p10 False-negative rate \eqn{p_{10} = \Pr(\hat Z = 0 \mid Z = 1)}
#'   (\code{K = 2} only). Auto-extracted from \code{Pi} when supplied.
#' @param pi_z Latent prevalence: a scalar
#'   \eqn{\pi = \Pr(Z = 1)} (\code{K = 2}) or a \eqn{K}-vector
#'   \eqn{(\pi_0, \dots, \pi_{K-1})} of class probabilities. If
#'   \code{NULL}, \code{pi_z} is estimated from \code{z_hat} and
#'   \code{Pi} via \eqn{\hat\pi = \Pi^{-1} \hat\pi_{\text{obs}}}
#'   (clamped to \eqn{[0.01, 0.99]}).
#' @param Pi The \eqn{K \times K} misclassification matrix
#'   \eqn{\Pi_{j\ell} = \Pr(\hat Z = j - 1 \mid Z = \ell - 1)}; columns
#'   must sum to 1. Extracted automatically from the \code{mc()} term
#'   when using the formula interface.
#' @param K Number of latent categories. Auto-detected from \code{Pi}
#'   or, failing that, from the levels of \code{z_hat}.
#' @param c1 Pre-computed misclassification constant
#'   \eqn{c_1 = p_{01}(1 - \pi)} (binary case, advanced use). Computed
#'   from \code{p01} and \code{pi_z} when \code{NULL}.
#' @param c2 Pre-computed misclassification constant
#'   \eqn{c_2 = p_{01}(1 - \pi) - p_{10} \pi} (binary case, advanced
#'   use). Computed from \code{p01}, \code{p10} and \code{pi_z} when
#'   \code{NULL}.
#' @param iterate Logical. If \code{TRUE}, iterate the BCA/BCM updates
#'   to convergence (the iterated BCM solves the corrected score
#'   equation; iterated BCA performs Newton steps on the additive
#'   correction).
#' @param jacobian Character. Either \code{"analytical"} (default) or
#'   \code{"numerical"}. Controls how the drift Jacobian
#'   \eqn{M = \partial m / \partial \psi} is evaluated for the
#'   multicategory (\eqn{K \ge 3}) BCM and CS estimators. The analytical
#'   form is exact and faster; the numerical form (forward differences)
#'   is provided as a fallback for cross-checking. The binary case
#'   (\eqn{K = 2}) always uses the analytical Jacobian.
#' @param fix_omega Logical. If \code{TRUE}, fix the mixture weights
#'   in the one-step estimator at the values implied by the supplied
#'   misclassification parameters; if \code{FALSE} (default), they are
#'   estimated jointly with \eqn{\psi} via a softmax parameterization.
#' @param vcov_corrected Logical. If \code{TRUE}, the BCA/BCM variance
#'   accounts for the additional uncertainty due to plug-in estimation
#'   of \eqn{\hat m(\hat\psi)} (joint score-and-drift sandwich); if
#'   \code{FALSE} (default), the drifting-regime asymptotic variance
#'   \eqn{A^{-1} C A^{-1}} is used. Both share the same point estimate.
#' @param weights Optional positive frequency weights of length \eqn{n}.
#' @param J Number of response categories (multinomial family only;
#'   auto-detected).
#' @param homoskedastic Logical. For one-step Gaussian fits, assume a
#'   single residual variance across mixture components; if
#'   \code{FALSE}, separate variances are estimated.
#' @param optim_control List of control parameters passed to
#'   \code{\link[stats]{nlminb}} for the one-step estimator.
#' @param z_hat Integer vector of observed proxy values
#'   \eqn{\hat Z_i \in \{0, \dots, K-1\}} (matrix interface only).
#' @param x Covariate matrix \eqn{[x_1, \dots, x_n]^\top} of dimension
#'   \eqn{n \times r}; should include an intercept column when one is
#'   desired (matrix interface only).
#'
#' @return An object of class \code{"mcglm"} with components:
#'   \describe{
#'     \item{coefficients}{Named list of coefficient vectors
#'       \eqn{\hat\psi}, one per method.}
#'     \item{vcov}{Named list of asymptotic variance-covariance
#'       matrices, one per method (sandwich estimators from the
#'       paper's theorems).}
#'     \item{se}{Named list of standard errors per method.}
#'     \item{naive_fit}{The \code{\link[stats]{glm}} object from the
#'       naive fit (\code{NULL} for multinomial).}
#'     \item{method, family, K, n, p}{Metadata.}
#'     \item{weights}{Frequency weights used (\code{NULL} if
#'       unweighted).}
#'     \item{loglik_onestep, vcov_onestep}{One-step log-likelihood
#'       and variance (when \code{onestep} is fit).}
#'     \item{call, formula}{The matched call and the formula
#'       (\code{NULL} if the matrix interface was used).}
#'   }
#'   See \code{\link{summary.mcglm}}, \code{\link{vcov.mcglm}},
#'   \code{\link{confint.mcglm}}, \code{\link{predict.mcglm}}, and
#'   \code{\link{logLik.mcglm}} for downstream methods.
#'
#' @references
#' Battaglia, L., Christensen, T., Hansen, S. and Sacher, S. (2025).
#' Inference for regression with variables generated by AI or machine
#' learning. \emph{arXiv preprint arXiv:2402.15585}.
#'
#' @seealso \code{\link{simex}} for SIMEX and MC-SIMEX corrections;
#'   \code{\link{mc}} for the formula-interface marker;
#'   \code{\link{summary.mcglm}}, \code{\link{vcov.mcglm}},
#'   \code{\link{confint.mcglm}}.
#'
#' @examples
#' # --- Formula interface (binary misclassification, Poisson) ---
#' set.seed(42)
#' n  <- 2000
#' Pi <- matrix(c(0.90, 0.10,    # Pi[2,1] = p01 = 0.10
#'                0.15, 0.85),   # Pi[1,2] = p10 = 0.15
#'              nrow = 2, byrow = FALSE)
#' z      <- rbinom(n, 1, 0.4)               # latent Z, prevalence pi = 0.4
#' z_hat  <- z
#' z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)        # flip 0 -> 1
#' z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)    # flip 1 -> 0
#' x1 <- rnorm(n)
#' y  <- rpois(n, exp(0.8 * z + -0.5 + 0.7 * x1))
#' df <- data.frame(y = y, z = factor(z_hat), x1 = x1)
#'
#' fit <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
#'              method = c("naive", "bca", "bcm", "cs"))
#' fit
#' summary(fit)
#' confint(fit, method = "cs")
#'
#' # --- Matrix interface, supplying p01/p10/pi_z directly ---
#' x_mat <- cbind(1, x1)
#' fit2 <- mcglm(y, z_hat = as.integer(z_hat), x = x_mat,
#'               family = "poisson", method = c("naive", "cs"),
#'               p01 = 0.10, p10 = 0.15, pi_z = 0.4)
#'
#' # --- Multicategory misclassification (K = 3) ---
#' \donttest{
#' set.seed(3)
#' n   <- 1000
#' Pi3 <- matrix(c(0.8, 0.1, 0.1,
#'                 0.1, 0.8, 0.1,
#'                 0.1, 0.1, 0.8), 3, 3, byrow = TRUE)
#' pi3 <- c(0.5, 0.3, 0.2)
#' Z   <- sample(0:2, n, replace = TRUE, prob = pi3)
#' Zh  <- vapply(Z, function(z) sample(0:2, 1, prob = Pi3[, z + 1]), 0L)
#' x1  <- rnorm(n)
#' y   <- rpois(n, exp(c(0, 1.0, -0.9)[Z + 1] + 0.8 - 0.7 * x1))
#' fit3 <- mcglm(y, z_hat = Zh, x = cbind(1, x1), family = "poisson",
#'               method = c("naive", "cs"), Pi = Pi3, pi_z = pi3)
#' coef(fit3, method = "cs")
#' }
#'
#' @export
mcglm <- function(formula, data = NULL, family = "poisson",
                  method = c("naive", "bca", "bcm", "cs"),
                  p01 = NULL, p10 = NULL, pi_z = NULL,
                  Pi = NULL, K = NULL,
                  c1 = NULL, c2 = NULL,
                  iterate = FALSE,
                  jacobian = c("analytical", "numerical"),
                  fix_omega = FALSE,
                  vcov_corrected = FALSE,
                  weights = NULL,
                  J = NULL,
                  homoskedastic = TRUE,
                  optim_control = list(),
                  z_hat = NULL, x = NULL) {

  cl <- match.call()
  jacobian <- match.arg(jacobian)
  formula_obj <- NULL

  # --- Dispatch: formula vs matrix interface ---
  if (inherits(formula, "formula")) {
    formula_obj <- formula
    caller_env <- parent.frame()
    parsed <- .mcglm_parse_formula(formula, data, caller_env)
    y     <- parsed$y
    z_hat <- parsed$z_hat
    x     <- parsed$x
    if (is.null(Pi) && !is.null(parsed$Pi)) Pi <- parsed$Pi
    if (is.null(K))  K  <- parsed$K
  } else {
    # Matrix interface: first arg is y
    y <- formula
  }

  # --- Derive p01/p10 from 2x2 Pi ---
  if (!is.null(Pi) && nrow(Pi) == 2L) {
    if (is.null(p01)) p01 <- Pi[2, 1]
    if (is.null(p10)) p10 <- Pi[1, 2]
  }

  # --- Estimate pi_z from observed data + Pi if not supplied ---
  if (is.null(pi_z) && !is.null(Pi) && !is.null(z_hat)) {
    pi_z <- .mcglm_estimate_pi_z(z_hat, Pi)
  }

  # --- Derive c1/c2 from (p01, p10, pi_z) if not supplied directly ---
  if (is.null(c1) && !is.null(p01) && !is.null(pi_z))
    c1 <- p01 * (1 - pi_z)
  if (is.null(c2) && !is.null(p01) && !is.null(p10) && !is.null(pi_z))
    c2 <- p01 * (1 - pi_z) - p10 * pi_z

  # --- Core fitting ---
  out <- .mcglm_fit(y = y, z_hat = z_hat, x = x, family = family,
                    method = method, p01 = p01, p10 = p10, pi_z = pi_z,
                    Pi = Pi, K = K, c1 = c1, c2 = c2, iterate = iterate,
                    jacobian = jacobian,
                    fix_omega = fix_omega,
                    vcov_corrected = vcov_corrected,
                    weights = weights, J = J,
                    homoskedastic = homoskedastic,
                    optim_control = optim_control)
  out$call    <- cl
  out$formula <- formula_obj
  out
}


#' Parse mcglm formula extracting mc() term
#' @keywords internal
.mcglm_parse_formula <- function(formula, data, env) {
  if (is.null(data))
    stop("'data' argument is required when using the formula interface.")

  lhs <- formula[[2]]
  rhs <- formula[[3]]

  # Extract response
  y_name <- deparse(lhs)
  y <- data[[y_name]]
  if (is.null(y))
    stop("Response variable '", y_name, "' not found in data.")


  # Walk RHS to find mc() term and other covariates
  mc_info <- NULL
  other_terms <- character(0)

  .walk_mcglm_rhs <- function(node) {
    if (is.call(node)) {
      op <- deparse(node[[1]])
      if (op == "+") {
        .walk_mcglm_rhs(node[[2]])
        .walk_mcglm_rhs(node[[3]])
        return(invisible(NULL))
      }
      if (op == "mc") {
        var_name <- deparse(node[[2]])
        if (length(node) >= 3L) {
          mat_val <- eval(node[[3]], data, env)
          mc_info <<- list(variable = var_name, Pi = as.matrix(mat_val))
        } else {
          mc_info <<- list(variable = var_name, Pi = NULL)
        }
        return(invisible(NULL))
      }
    }
    # Everything else is a regular term
    other_terms <<- c(other_terms, deparse(node))
    invisible(NULL)
  }

  .walk_mcglm_rhs(rhs)

  if (is.null(mc_info))
    stop("Formula must contain exactly one mc() term, e.g., y ~ mc(z, Pi) + x")

  # Get z_hat from data
  z_var <- data[[mc_info$variable]]
  if (is.null(z_var))
    stop("mc() variable '", mc_info$variable, "' not found in data.")

  Pi <- mc_info$Pi

  # Convert z to 0-based integer
  if (is.factor(z_var)) {
    z_hat <- as.integer(z_var) - 1L
  } else {
    z_hat <- as.integer(z_var)
  }

  # Determine K
  K <- if (!is.null(Pi)) nrow(Pi) else length(unique(z_hat))

  # Validate Pi when provided
  if (!is.null(Pi)) {
    if (ncol(Pi) != K)
      stop("Pi must be a square K x K matrix (got ", nrow(Pi), " x ", ncol(Pi), ").")
    cs <- colSums(Pi)
    if (any(abs(cs - 1) > 1e-6))
      stop("Columns of Pi must sum to 1 (got: ",
           paste(round(cs, 4), collapse = ", "), ").")
  }

  # Build x matrix from remaining terms
  if (length(other_terms) == 0) {
    # No other covariates: just intercept
    x <- matrix(1, nrow = nrow(data), ncol = 1)
  } else {
    rhs_formula <- as.formula(paste("~", paste(other_terms, collapse = " + ")),
                              env = env)
    x <- model.matrix(rhs_formula, data = data)
  }

  list(y = as.numeric(y), z_hat = z_hat, x = x, Pi = Pi, K = K)
}


#' Core mcglm fitting (internal, called by mcglm after dispatch)
#' @keywords internal
.mcglm_fit <- function(y, z_hat, x, family = "poisson",
                       method = c("naive", "bca", "bcm", "cs"),
                       p01 = NULL, p10 = NULL, pi_z = NULL,
                       Pi = NULL, K = NULL,
                       c1 = NULL, c2 = NULL,
                       iterate = FALSE,
                       jacobian = c("analytical", "numerical"),
                       fix_omega = FALSE,
                       vcov_corrected = FALSE,
                       weights = NULL,
                       J = NULL,
                       homoskedastic = TRUE,
                       optim_control = list()) {

  # --- input validation ---
  y     <- as.numeric(y)
  z_hat <- as.integer(z_hat)
  x     <- as.matrix(x)
  n     <- length(y)
  jacobian <- match.arg(jacobian)
  method <- match.arg(method, c("naive", "bca", "bcm", "cs", "onestep"),
                      several.ok = TRUE)

  stopifnot(nrow(x) == n, length(z_hat) == n)

  # --- validate frequency weights ---
  wt <- NULL
  if (!is.null(weights)) {
    wt <- as.numeric(weights)
    if (length(wt) != n)
      stop("weights must have length n (", n, "), got ", length(wt))
    if (any(wt <= 0))
      stop("weights must be positive")
  }

  # --- multinomial family ---
  is_multinomial <- is.character(family) && family == "multinomial"
  if (is_multinomial) {
    unsupported <- setdiff(method, c("naive", "onestep"))
    if (length(unsupported) > 0)
      stop("For multinomial family, only 'naive' and 'onestep' methods are ",
           "supported. Unsupported: ", paste(unsupported, collapse = ", "))
    y <- as.integer(y)
    if (is.null(J)) J <- length(unique(y))
    if (any(y < 0L) || any(y >= J))
      stop("For multinomial, y must be integers in {0, ..., J-1}")
  }

  # Determine binary vs multicategory
  if (is.null(K)) {
    unique_z <- sort(unique(z_hat))
    if (is.null(Pi) && all(unique_z %in% c(0L, 1L))) {
      K <- 2L
    } else {
      K <- length(unique_z)
    }
  }
  is_binary <- (K == 2L)

  # --- validate misclassification parameters ---
  needs_correction <- any(method %in% c("bca", "bcm", "cs", "onestep"))
  if (needs_correction) {
    if (is_binary) {
      needs_bca_bcm_cs <- any(method %in% c("bca", "bcm", "cs"))
      needs_onestep    <- "onestep" %in% method
      if (needs_bca_bcm_cs && (is.null(c1) || is.null(c2)))
        stop("For binary bca/bcm/cs, supply c1 and c2 (or Pi, or p01/p10/pi_z).")
      if (needs_onestep && fix_omega && (is.null(p01) || is.null(p10) || is.null(pi_z)))
        stop("For onestep with fix_omega=TRUE, supply p01, p10, and pi_z.")
    } else {
      needs_bca_bcm_cs <- any(method %in% c("bca", "bcm", "cs"))
      needs_onestep    <- "onestep" %in% method
      if (needs_bca_bcm_cs) {
        if (is.null(Pi) || is.null(pi_z))
          stop("For multicategory bca/bcm/cs, supply Pi and pi_z.")
        Pi   <- as.matrix(Pi)
        pi_z <- as.numeric(pi_z)
        stopifnot(nrow(Pi) == K, ncol(Pi) == K, length(pi_z) == K)
      }
      if (needs_onestep && fix_omega && (is.null(Pi) || is.null(pi_z)))
        stop("For multicategory onestep with fix_omega=TRUE, supply Pi and pi_z.")
    }
  }

  # --- build xi_hat once ---
  if (!is_multinomial) {
    xi_hat <- .mcglm_build_xi_hat(z_hat, x, K)
  }

  # --- fit ---
  results <- list()
  onestep_vcov   <- NULL
  onestep_loglik <- NULL
  naive_glm_fit  <- NULL

  if (is_multinomial) {
    naive_mn <- .mcglm_fit_naive_multinomial(y, z_hat, x, J, K, wt = wt)
    psi      <- naive_mn$coefficients
    results$naive <- psi
    p_total  <- length(psi)

    if ("onestep" %in% method) {
      os <- .mcglm_fit_onestep_multinomial(y, z_hat, x, J, K,
                                           p01 = p01, p10 = p10, pi_z = pi_z,
                                           Pi = Pi, fix_omega = fix_omega,
                                           optim_control = optim_control,
                                           wt = wt)
      results$onestep <- os$coefficients
      onestep_vcov    <- os$vcov
      onestep_loglik  <- os$loglik
    }

  } else if (is_binary) {
    naive <- .mcglm_fit_naive_bin(y, xi_hat, family, wt = wt)
    psi   <- naive$coefficients
    results$naive <- psi
    p_total <- length(psi)
    naive_glm_fit <- naive$glm_fit

    if ("bca" %in% method)
      results$bca <- .mcglm_fit_bca_bin(psi, y, xi_hat, x, family, c1, c2,
                                         iterate = iterate, wt = wt)
    if ("bcm" %in% method)
      results$bcm <- .mcglm_fit_bcm_bin(psi, y, xi_hat, x, family, c1, c2,
                                         iterate = iterate, wt = wt)
    if ("cs"  %in% method)
      results$cs  <- .mcglm_fit_cs_bin(psi, y, xi_hat, x, family, c1, c2,
                                        wt = wt)

    if ("onestep" %in% method) {
      os <- .mcglm_fit_onestep_bin(y, z_hat, x, family, p01, p10, pi_z,
                                   fix_omega = fix_omega,
                                   homoskedastic = homoskedastic,
                                   optim_control = optim_control, wt = wt)
      results$onestep <- os$coefficients
      onestep_vcov    <- os$vcov
      onestep_loglik  <- os$loglik
    }

  } else {
    naive <- .mcglm_fit_naive_multi(y, xi_hat, family, wt = wt)
    psi   <- naive$coefficients
    results$naive <- psi
    p_total <- length(psi)
    naive_glm_fit <- naive$glm_fit

    if ("bca" %in% method)
      results$bca <- .mcglm_fit_bca_multi(psi, y, xi_hat, z_hat, x, K, family,
                                           Pi, pi_z, iterate = iterate, wt = wt)
    if ("bcm" %in% method)
      results$bcm <- .mcglm_fit_bcm_multi(psi, y, xi_hat, z_hat, x, K, family,
                                           Pi, pi_z, iterate = iterate, wt = wt,
                                           jacobian = jacobian)
    if ("cs"  %in% method)
      results$cs  <- .mcglm_fit_cs_multi(psi, y, xi_hat, z_hat, x, K, family,
                                          Pi, pi_z, wt = wt,
                                          jacobian = jacobian)

    if ("onestep" %in% method) {
      os <- .mcglm_fit_onestep_multi(y, z_hat, x, K, family, Pi, pi_z,
                                     fix_omega = fix_omega,
                                     homoskedastic = homoskedastic,
                                     optim_control = optim_control, wt = wt)
      results$onestep <- os$coefficients
      onestep_vcov    <- os$vcov
      onestep_loglik  <- os$loglik
    }
  }

  # --- parameter names ---
  if (is_multinomial) {
    s <- K - 1
    r <- ncol(x)
    block_size <- s + r
    nms <- character(p_total)
    for (jj in seq_len(J - 1)) {
      offset <- (jj - 1) * block_size
      if (s == 1) {
        nms[offset + 1] <- paste0("y", jj, ":gamma")
      } else if (s > 1) {
        nms[offset + seq_len(s)] <- paste0("y", jj, ":gamma", seq_len(s))
      }
      nms[offset + s + seq_len(r)] <- paste0("y", jj, ":alpha", seq_len(r) - 1)
    }
  } else if (is_binary) {
    nms <- c("gamma", paste0("alpha", seq_len(ncol(x)) - 1))
  } else {
    nms <- c(paste0("gamma", seq_len(K - 1)),
             paste0("alpha", seq_len(ncol(x)) - 1))
  }
  results <- lapply(results, function(v) {
    v <- as.numeric(v)
    names(v) <- nms
    v
  })

  # --- per-method asymptotic variance ---
  vcov_list <- list()
  if (!is_multinomial) {
    for (nm in names(results)) {
      psi_nm <- unname(results[[nm]])
      V <- tryCatch(
        if (nm == "naive") {
          if (is_binary)
            .mcglm_vcov_naive(psi_nm, y, xi_hat, family, wt = wt)
          else
            .mcglm_vcov_naive_multi(psi_nm, y, xi_hat, z_hat, x, K, family,
                                    wt = wt)
        } else if (nm %in% c("bca", "bcm")) {
          have_probs <- (!is.null(p01) && !is.null(p10) && !is.null(pi_z)) ||
                        (!is.null(c1) && !is.null(c2))
          if (is_binary) {
            .mcglm_vcov_bc_bin(psi_nm, y, xi_hat, x, family,
                               p01 = p01, p10 = p10, pi_z = pi_z,
                               c1 = c1, c2 = c2,
                               psi_naive = unname(results$naive),
                               type = nm,
                               corrected = vcov_corrected && have_probs,
                               wt = wt)
          } else {
            have_multi <- !is.null(Pi) && !is.null(pi_z)
            .mcglm_vcov_bc_multi(psi_nm, y, xi_hat, z_hat, x, K, family,
                                 Pi = Pi, pi_z = pi_z,
                                 psi_naive = unname(results$naive),
                                 type = nm,
                                 corrected = vcov_corrected && have_multi,
                                 wt = wt,
                                 jacobian = jacobian)
          }
        } else if (nm == "cs") {
          if (is_binary)
            .mcglm_vcov_cs_bin(psi_nm, y, xi_hat, x, family, p01, p10, pi_z,
                               c1 = c1, c2 = c2,
                               wt = wt)
          else
            .mcglm_vcov_cs_multi(psi_nm, y, xi_hat, z_hat, x, K, family,
                                 Pi, pi_z, wt = wt,
                                 jacobian = jacobian)
        } else if (nm == "onestep") {
          onestep_vcov
        },
        error = function(e) {
          warning("Variance computation failed for method '", nm, "': ",
                  conditionMessage(e))
          NULL
        }
      )
      if (!is.null(V)) {
        dimnames(V) <- list(nms, nms)
      }
      vcov_list[[nm]] <- V
    }
  } else if (!is.null(onestep_vcov)) {
    vcov_list$onestep <- onestep_vcov
  }

  se_list <- lapply(vcov_list, function(V) {
    if (is.null(V)) return(NULL)
    s <- sqrt(pmax(diag(V), 0))
    names(s) <- rownames(V)
    s
  })

  out <- list(
    coefficients = results,
    vcov         = vcov_list,
    se           = se_list,
    naive_fit    = naive_glm_fit,
    method       = method,
    family       = family,
    K            = K,
    n            = n,
    p            = p_total,
    weights      = wt,
    y            = y,
    z_hat        = z_hat,
    x            = x,
    is_multinomial = is_multinomial
  )
  if (!is_multinomial) out$xi_hat <- xi_hat
  if (is_multinomial)  out$J     <- J

  if (!is.null(onestep_vcov)) {
    out$vcov_onestep   <- onestep_vcov  # kept for backwards compatibility
    out$loglik_onestep <- onestep_loglik
  }

  structure(out, class = "mcglm")
}


#' Default method for an \code{mcglm} object
#' @keywords internal
.mcglm_default_method <- function(object) {
  if (!is.null(object$method) && length(object$method))
    return(object$method[length(object$method)])
  names(object$coefficients)[1]
}

.mcglm_family_string <- function(object) {
  if (is.character(object$family)) return(object$family)
  if (inherits(object$family, "family")) return(object$family$family)
  as.character(object$family)
}


#' Print method for \code{mcglm} fits
#'
#' Mirrors the layout of \code{print.glm}: shows the call, the family,
#' a coefficient table (one column per estimation method), the residual
#' degrees of freedom and (when available) the deviance and AIC of the
#' naive GLM fit.
#' @param x An object of class \code{"mcglm"}.
#' @param digits Number of digits to display.
#' @param ... Unused.
#' @export
print.mcglm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n")
  if (!is.null(x$call)) print(x$call) else cat("mcglm(...)\n")
  cat("\nFamily: ", .mcglm_family_string(x),
      "  |  n = ", x$n, ", K = ", x$K,
      ", p = ", x$p, "\n", sep = "")
  cat("Methods: ", paste(x$method, collapse = ", "), "\n\n", sep = "")

  cat("Coefficients:\n")
  coefs <- do.call(cbind, x$coefficients)
  colnames(coefs) <- toupper(names(x$coefficients))
  print.default(format(coefs, digits = digits), print.gap = 2L,
                quote = FALSE)

  if (!is.null(x$naive_fit)) {
    rdf <- x$naive_fit$df.residual
    cat("\nDegrees of Freedom:", x$naive_fit$df.null, "Total (i.e. Null); ",
        rdf, "Residual\n")
    cat("Null Deviance:    ",
        format(signif(x$naive_fit$null.deviance, digits)), "\n")
    cat("Residual Deviance:",
        format(signif(x$naive_fit$deviance, digits)),
        " | AIC (naive):", format(signif(stats::AIC(x$naive_fit), digits)),
        "\n")
  }
  if (!is.null(x$loglik_onestep)) {
    cat("logLik (onestep): ", format(signif(x$loglik_onestep, digits)), "\n")
  }
  invisible(x)
}


#' Extract coefficients from an \code{mcglm} fit
#'
#' @param object An \code{mcglm} object.
#' @param method Optional character: which estimator to extract
#'   (\code{"naive"}, \code{"bca"}, \code{"bcm"}, \code{"cs"},
#'   \code{"onestep"}). If \code{NULL}, returns a named list of all
#'   coefficient vectors.
#' @param ... Unused.
#' @export
coef.mcglm <- function(object, method = NULL, ...) {
  if (is.null(method)) return(object$coefficients)
  if (!method %in% names(object$coefficients))
    stop("Method '", method, "' was not fit. Available: ",
         paste(names(object$coefficients), collapse = ", "))
  object$coefficients[[method]]
}


#' Variance-covariance matrix of an \code{mcglm} fit
#'
#' @param object An \code{mcglm} object.
#' @param method Optional character: which estimator to use. Defaults to the
#'   last method fit (e.g. \code{"onestep"} when present, else the most
#'   refined correction available).
#' @param ... Unused.
#' @return A \eqn{p \times p} sandwich variance-covariance matrix.
#' @export
vcov.mcglm <- function(object, method = NULL, ...) {
  if (is.null(object$vcov)) {
    # backwards-compat (older fits stored only vcov_onestep)
    return(object$vcov_onestep)
  }
  if (is.null(method)) method <- .mcglm_default_method(object)
  if (!method %in% names(object$vcov))
    stop("No variance available for method '", method, "'.")
  object$vcov[[method]]
}


#' Summary of an \code{mcglm} fit
#'
#' Returns a Wald-style coefficient table (estimate, standard error,
#' z-value, p-value) for each fitted method, using the asymptotic
#' sandwich variance derived in Battaglia et al. (2025).
#' @param object An \code{mcglm} object.
#' @param ... Unused.
#' @export
summary.mcglm <- function(object, ...) {
  fam_str <- .mcglm_family_string(object)
  ans <- list(
    call    = object$call,
    family  = fam_str,
    n       = object$n,
    p       = object$p,
    K       = object$K,
    method  = object$method,
    coefficients = list(),
    deviance     = if (!is.null(object$naive_fit)) object$naive_fit$deviance,
    df.residual  = if (!is.null(object$naive_fit)) object$naive_fit$df.residual,
    aic_naive    = if (!is.null(object$naive_fit)) stats::AIC(object$naive_fit),
    loglik_onestep = object$loglik_onestep
  )
  for (nm in names(object$coefficients)) {
    est <- object$coefficients[[nm]]
    V   <- if (!is.null(object$vcov)) object$vcov[[nm]] else NULL
    if (!is.null(V)) {
      se <- sqrt(pmax(diag(V), 0))
      z  <- est / se
      pv <- 2 * stats::pnorm(-abs(z))
      tab <- cbind(Estimate = est,
                   `Std. Error` = se,
                   `z value` = z,
                   `Pr(>|z|)` = pv)
    } else {
      tab <- cbind(Estimate = est)
    }
    ans$coefficients[[nm]] <- tab
  }
  class(ans) <- "summary.mcglm"
  ans
}


#' @export
print.summary.mcglm <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  cat("\nCall:\n")
  if (!is.null(x$call)) print(x$call) else cat("mcglm(...)\n")
  cat("\nFamily: ", x$family,
      "  |  n = ", x$n, ", K = ", x$K, ", p = ", x$p, "\n", sep = "")
  cat("Methods: ", paste(x$method, collapse = ", "), "\n", sep = "")

  for (nm in x$method) {
    cat("\n--- ", toupper(nm), " ---\n", sep = "")
    tab <- x$coefficients[[nm]]
    if (is.null(tab)) next
    has_p <- ncol(tab) >= 4L
    stats::printCoefmat(tab, digits = digits,
                        has.Pvalue = has_p, P.values = has_p,
                        signif.stars = signif.stars && has_p,
                        cs.ind = if (has_p) 1:2 else 1L,
                        tst.ind = if (has_p) 3L else integer(0),
                        na.print = "NA")
  }

  if (!is.null(x$deviance)) {
    cat("\nResidual deviance (naive): ",
        format(signif(x$deviance, digits)),
        "  on ", x$df.residual, " degrees of freedom\n", sep = "")
    cat("AIC (naive): ", format(signif(x$aic_naive, digits)), "\n", sep = "")
  }
  if (!is.null(x$loglik_onestep)) {
    cat("logLik (onestep): ", format(signif(x$loglik_onestep, digits)), "\n",
        sep = "")
  }
  if (length(x$coefficients) > 1L) {
    cat("\nBias correction (difference from naive):\n")
    naive <- x$coefficients[["naive"]][, "Estimate"]
    others <- setdiff(names(x$coefficients), "naive")
    diffs <- sapply(others, function(o)
      x$coefficients[[o]][, "Estimate"] - naive)
    if (!is.matrix(diffs)) diffs <- matrix(diffs, nrow = 1L,
                                           dimnames = list(names(naive), others))
    print.default(format(diffs, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  invisible(x)
}


#' Standard errors of coefficient estimates
#'
#' @param object An \code{mcglm} object.
#' @param method Optional method name.
#' @return Named numeric vector of standard errors.
#' @export
se.mcglm <- function(object, method = NULL) {
  if (is.null(object$se)) return(sqrt(diag(vcov(object, method = method))))
  if (is.null(method)) method <- .mcglm_default_method(object)
  object$se[[method]]
}


#' Wald confidence intervals for an \code{mcglm} fit
#'
#' @param object An \code{mcglm} object.
#' @param parm Optional names or indices of parameters to return.
#' @param level Coverage level.
#' @param method Optional estimation method (defaults to the last one fit,
#'   e.g. \code{"onestep"} if present).
#' @param ... Unused.
#' @export
confint.mcglm <- function(object, parm, level = 0.95, method = NULL, ...) {
  if (is.null(method)) method <- .mcglm_default_method(object)
  est <- coef(object, method = method)
  V   <- vcov(object, method = method)
  if (is.null(V))
    stop("No variance available for method '", method, "'.")
  se   <- sqrt(pmax(diag(V), 0))
  alpha <- (1 - level) / 2
  q    <- stats::qnorm(1 - alpha)
  ci   <- cbind(est - q * se, est + q * se)
  pct  <- paste0(format(c(alpha, 1 - alpha) * 100,
                        trim = TRUE, scientific = FALSE,
                        digits = 3), " %")
  colnames(ci) <- pct
  rownames(ci) <- names(est)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}


#' Fitted values for an \code{mcglm} fit
#'
#' @param object An \code{mcglm} object.
#' @param method Estimation method.
#' @param ... Unused.
#' @export
fitted.mcglm <- function(object, method = NULL, ...) {
  if (isTRUE(object$is_multinomial))
    stop("fitted() not yet implemented for multinomial mcglm fits.")
  if (is.null(method)) method <- .mcglm_default_method(object)
  psi <- coef(object, method = method)
  fam <- .mcglm_get_link_funs(object$family)
  eta <- as.numeric(object$xi_hat %*% unname(psi))
  setNames(fam$mu(eta), seq_along(eta))
}


#' Linear-predictor or response predictions
#'
#' Currently supports in-sample prediction only.
#' @param object An \code{mcglm} object.
#' @param newdata Currently unused. Predictions are computed at the in-sample
#'   proxy design \eqn{\hat\xi_i}.
#' @param type One of \code{"link"}, \code{"response"}.
#' @param method Estimation method.
#' @param ... Unused.
#' @export
predict.mcglm <- function(object, newdata = NULL,
                          type = c("link", "response"),
                          method = NULL, ...) {
  type <- match.arg(type)
  if (!is.null(newdata))
    stop("predict() with 'newdata' not yet supported for mcglm fits.")
  if (isTRUE(object$is_multinomial))
    stop("predict() not yet implemented for multinomial mcglm fits.")
  if (is.null(method)) method <- .mcglm_default_method(object)
  psi <- coef(object, method = method)
  eta <- as.numeric(object$xi_hat %*% unname(psi))
  if (type == "link") return(eta)
  fam <- .mcglm_get_link_funs(object$family)
  fam$mu(eta)
}


#' Residuals for an \code{mcglm} fit
#'
#' @param object An \code{mcglm} object.
#' @param method Estimation method.
#' @param type One of \code{"response"}, \code{"pearson"}, \code{"deviance"}.
#' @param ... Unused.
#' @export
residuals.mcglm <- function(object, method = NULL,
                            type = c("response", "pearson", "deviance"),
                            ...) {
  type <- match.arg(type)
  if (isTRUE(object$is_multinomial))
    stop("residuals() not yet implemented for multinomial mcglm fits.")
  if (is.null(method)) method <- .mcglm_default_method(object)
  fam_funs <- .mcglm_get_link_funs(object$family)
  fam      <- fam_funs$family
  mu_val   <- fitted(object, method = method)
  r        <- object$y - mu_val
  if (type == "response") return(r)
  if (type == "pearson") return(r / sqrt(fam$variance(mu_val)))
  wt <- if (is.null(object$weights)) rep.int(1, length(r)) else object$weights
  d  <- fam$dev.resids(object$y, mu_val, wt = wt)
  sign(r) * sqrt(d)
}


#' Number of observations
#' @param object An \code{mcglm} object.
#' @param ... Unused.
#' @export
nobs.mcglm <- function(object, ...) object$n


#' Family of an \code{mcglm} fit
#' @param object An \code{mcglm} object.
#' @param ... Unused.
#' @export
family.mcglm <- function(object, ...) {
  .mcglm_get_link_funs(object$family)$family
}


#' Formula of an \code{mcglm} fit (NULL for matrix interface)
#' @param x An \code{mcglm} object.
#' @param ... Unused.
#' @export
formula.mcglm <- function(x, ...) x$formula


#' Proxy design matrix used by \code{mcglm}
#'
#' Returns the assembled \eqn{\hat\xi_i = (\hat d_i, x_i)} design matrix
#' (proxy-encoded misclassified covariate plus other covariates).
#' @param object An \code{mcglm} object.
#' @param ... Unused.
#' @export
model.matrix.mcglm <- function(object, ...) object$xi_hat


#' Log-likelihood of an \code{mcglm} fit
#'
#' Returns the naive-GLM log-likelihood for \code{method = "naive"} (the
#' only estimator with a closed-form likelihood for proxy data) or the
#' integrated mixture log-likelihood for \code{method = "onestep"}. Other
#' methods do not correspond to a likelihood; an error is thrown.
#' @param object An \code{mcglm} object.
#' @param method Estimation method.
#' @param ... Unused.
#' @export
logLik.mcglm <- function(object, method = NULL, ...) {
  if (is.null(method))
    method <- if ("onestep" %in% object$method) "onestep" else "naive"
  if (method == "onestep" && !is.null(object$loglik_onestep)) {
    val <- object$loglik_onestep
    df  <- length(coef(object, method = "onestep"))
    attr(val, "df")   <- df
    attr(val, "nobs") <- object$n
    class(val) <- "logLik"
    return(val)
  }
  if (method == "naive" && !is.null(object$naive_fit)) {
    return(stats::logLik(object$naive_fit))
  }
  stop("logLik not defined for method '", method, "'. ",
       "Available: ", paste(intersect(c("naive", "onestep"), object$method),
                            collapse = ", "))
}


#' AIC of an \code{mcglm} fit
#' @param object An \code{mcglm} object.
#' @param ... Additional arguments (ignored).
#' @param k Penalty per parameter (default 2).
#' @export
AIC.mcglm <- function(object, ..., k = 2) {
  ll <- logLik(object)
  -2 * as.numeric(ll) + k * attr(ll, "df")
}


#' BIC of an \code{mcglm} fit
#' @param object An \code{mcglm} object.
#' @param ... Additional arguments (ignored).
#' @export
BIC.mcglm <- function(object, ...) {
  ll <- logLik(object)
  -2 * as.numeric(ll) + log(attr(ll, "nobs")) * attr(ll, "df")
}
