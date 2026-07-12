# ---------------------------------------------------------------------------
# nonprob_mc() -- mass-imputation estimator with a misclassification-corrected
# outcome model (mismeasured::mcglm), returning a nonprob-class object so the
# nonprobsvy S3 methods work unchanged.
# ---------------------------------------------------------------------------

#' Mass imputation with a misclassified covariate in the outcome model
#'
#' Estimates the population mean of an outcome observed only in a
#' non-probability sample by mass imputation onto a probability sample,
#' correcting the outcome model for a misclassified discrete covariate.
#' The interface mirrors \code{\link[nonprobsvy]{nonprob}} (mass-imputation
#' arguments only); the misclassified covariate is declared with the
#' \code{\link[mismeasured]{mc}} formula term, e.g.
#' \code{y ~ mc(z, Pi) + x1 + x2}.
#'
#' The outcome model is fit by \code{\link[mismeasured]{mcglm}} on the
#' non-probability sample, where only the misclassified proxy of \code{z}
#' is available. Predictions are made on the probability sample, which is
#' assumed to carry the \emph{true} \code{z} (see \sQuote{Details}), and
#' the population mean is the survey-weighted mean of these predictions.
#'
#' @details
#' \strong{Estimator.} With \eqn{\hat\psi} from the outcome model chosen by
#' \code{mc_method} and \eqn{w_i} the design weights,
#' \deqn{\hat{\bar Y} = \frac{\sum_{i \in S_P} w_i \,
#'   \mu(\hat\psi^\top \xi_i)}{\sum_{i \in S_P} w_i},}
#' computed via \code{\link[survey]{svymean}} (Hajek form, matching
#' \code{nonprob()}'s mass-imputation estimator).
#'
#' \strong{Variance.} \code{var_prob} is the design-based variance of the
#' weighted mean of predictions, from \pkg{survey} (any design supported
#' by \code{svydesign()}). \code{var_nonprob} accounts for the estimation
#' of \eqn{\hat\psi} via the corrected-score linearization
#' \deqn{\widehat{V}_{NP} = G_P^\top B \hat S B^\top G_P / n_{NP},}
#' with \eqn{B} and \eqn{\hat S} the \code{\link[sandwich]{bread}} and
#' \code{\link[sandwich]{meat}} of the \code{mcglm} fit and
#' \eqn{G_P = \hat N^{-1} \sum_{i \in S_P} w_i\, \dot\mu(\hat\psi^\top
#' \xi_i)\, \xi_i}. For \code{mc_method = "cs"} this is the variance from
#' the corrected-score mass-imputation theory; for \code{"naive"},
#' \code{"bca"} and \code{"bcm"} it uses the naive sandwich, their
#' asymptotic variance in the drifting regime.
#'
#' \strong{The probability sample must carry the true covariate.} The
#' column named inside \code{mc()} is taken at face value in
#' \code{svydesign$variables}: pass the true category there. If the
#' probability sample also only has a proxy, the returned estimator does
#' not correct for that (a marginalization extension is planned).
#'
#' \strong{Misclassification parameters.} Supplied exactly as in
#' \code{\link[mismeasured]{mcglm}}: a \eqn{K \times K} column-stochastic
#' matrix \code{Pi} (inside \code{mc()} or as an argument), or
#' \code{p01}/\code{p10}/\code{pi_z} in the binary case. They are treated
#' as known; the reported variance conditions on them.
#'
#' @param data A \code{data.frame} with the non-probability sample
#'   (outcome, proxy covariate, auxiliary variables).
#' @param outcome Formula for the outcome model with exactly one
#'   \code{\link[mismeasured]{mc}} term, e.g. \code{y ~ mc(z, Pi) + x1}.
#' @param svydesign A \code{\link[survey]{svydesign}} object describing the
#'   probability sample; its variables must include the covariates of
#'   \code{outcome} with the \code{mc()} variable holding the \emph{true}
#'   category.
#' @param pop_size Optional fixed population size; defaults to the sum of
#'   the design weights.
#' @param method_outcome Outcome-model method; only \code{"mcglm"} is
#'   available.
#' @param family_outcome Outcome family: \code{"gaussian"},
#'   \code{"binomial"} or \code{"poisson"} (canonical links, as required
#'   by the correction theory).
#' @param mc_method Which misclassification correction the estimator uses:
#'   \code{"cs"} (corrected score, default), \code{"bca"}, \code{"bcm"}
#'   or \code{"naive"}. All requested corrections are kept in the fitted
#'   outcome model for diagnostics.
#' @param p01,p10,pi_z,Pi Misclassification parameters, passed to
#'   \code{\link[mismeasured]{mcglm}}; see \sQuote{Details}.
#' @param case_weights Optional frequency weights for the non-probability
#'   sample.
#' @param control_outcome,control_inference Control lists, as in
#'   \code{\link[nonprobsvy]{nonprob}} (stored on the object; used by the
#'   \code{nonprobsvy} methods).
#' @param verbose Logical, print progress information.
#' @param se Logical, compute the variance (default \code{TRUE}).
#'
#' @return An object of classes \code{nonprob} and \code{nonprob_mc} with
#'   the same structure as a mass-imputation fit from
#'   \code{\link[nonprobsvy]{nonprob}}, so \code{print()},
#'   \code{summary()}, \code{confint()}, \code{nobs()}, \code{weights()}
#'   and \code{pop_size()} from \pkg{nonprobsvy} work unchanged. The
#'   \code{outcome} element stores the full \code{mcglm} fit
#'   (\code{fit$outcome[[1]]$model_fitted}), and \code{svydesign} is the
#'   updated design carrying the predictions in \code{.y_mc_hat} (useful
#'   for domain estimation with \code{\link[survey]{svyby}}; note the
#'   design-only variance there omits \code{var_nonprob}).
#'
#' @references
#' Kim, J. K., Park, S., Chen, Y. and Wu, C. (2021). Combining
#' non-probability and probability survey samples through mass imputation.
#' \emph{Journal of the Royal Statistical Society: Series A}, 184(3),
#' 941--963. \doi{10.1111/rssa.12696}
#'
#' @examples
#' \donttest{
#' library(survey)
#' set.seed(1)
#' N  <- 50000
#' x1 <- rnorm(N)
#' z  <- rbinom(N, 1, 0.4)
#' y  <- rpois(N, exp(0.8 * z - 0.5 + 0.7 * x1))
#' Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
#' z_hat <- ifelse(z == 1, 1 - rbinom(N, 1, 0.15), rbinom(N, 1, 0.10))
#'
#' # non-probability sample (proxy z), probability sample (true z, no y)
#' s_np  <- rbinom(N, 1, plogis(-2.5 + 0.5 * x1)) == 1
#' df_np <- data.frame(y = y, z = z_hat, x1 = x1)[s_np, ]
#' s_p   <- sample.int(N, 1000)
#' des   <- svydesign(ids = ~1, weights = rep(N / 1000, 1000),
#'                    data = data.frame(z = z[s_p], x1 = x1[s_p]))
#'
#' fit <- nonprob_mc(data = df_np, outcome = y ~ mc(z, Pi) + x1,
#'                   svydesign = des, family_outcome = "poisson")
#' fit
#' summary(fit)
#' cbind(naive = mean(y), corrected = fit$output$mean, true = mean(y))
#' }
#'
#' @seealso \code{\link[nonprobsvy]{nonprob}},
#'   \code{\link[mismeasured]{mcglm}}, \code{\link[mismeasured]{mc}}
#' @export
nonprob_mc <- function(data,
                       outcome = NULL,
                       svydesign = NULL,
                       pop_size = NULL,
                       method_outcome = c("mcglm"),
                       family_outcome = c("gaussian", "binomial", "poisson"),
                       mc_method = c("cs", "bca", "bcm", "naive"),
                       p01 = NULL,
                       p10 = NULL,
                       pi_z = NULL,
                       Pi = NULL,
                       case_weights = NULL,
                       control_outcome = control_out(),
                       control_inference = control_inf(),
                       verbose = FALSE,
                       se = TRUE) {

  cl <- match.call()
  method_outcome <- match.arg(method_outcome)
  family_outcome <- match.arg(family_outcome)
  mc_method      <- match.arg(mc_method)

  # ---- input checks --------------------------------------------------
  if (!inherits(outcome, "formula"))
    stop("'outcome' must be a formula with an mc() term, ",
         "e.g. y ~ mc(z, Pi) + x1.")
  if (!is.data.frame(data))
    stop("'data' must be a data.frame with the non-probability sample.")
  if (!inherits(svydesign, "survey.design"))
    stop("'svydesign' must be a survey::svydesign object ",
         "(the probability sample).")

  # parse the mc() outcome formula on both samples; the probability
  # sample carries no outcome, mc_parse_formula() allows that. Evaluate
  # the mc() matrix in the formula's own environment -- this function's
  # frame holds a (possibly NULL) 'Pi' argument that would shadow it.
  parts_nons <- mc_parse_formula(outcome, data, env = environment(outcome))
  if (is.null(Pi)) Pi <- parts_nons$Pi
  vars_rand  <- svydesign$variables
  mc_var     <- .mc_variable_name(outcome)
  if (!mc_var %in% names(vars_rand))
    stop("The mc() variable '", mc_var, "' was not found in ",
         "svydesign$variables. The probability sample must carry the ",
         "TRUE category of '", mc_var, "'.")
  parts_rand <- mc_parse_formula(outcome, vars_rand,
                                 env = environment(outcome),
                                 require_y = FALSE)
  if (any(parts_rand$z_hat < 0L) || any(parts_rand$z_hat >= parts_nons$K))
    stop("svydesign$variables$", mc_var, " has categories outside ",
         "0, ..., K-1 (K = ", parts_nons$K, ").")

  y_name <- all.vars(outcome[[2]])
  y_nons <- parts_nons$y
  n_nons <- length(y_nons)
  n_rand <- nrow(vars_rand)

  # ---- outcome model on the non-probability sample --------------------
  # fit the requested correction last so that it is mcglm's default
  # method for coef()/predict()/bread()/meat()
  methods_fit <- unique(c("naive", setdiff(c("bca", "bcm", "cs"), mc_method),
                          mc_method))
  if (verbose)
    message("Fitting mcglm (", family_outcome, ") with methods: ",
            paste(methods_fit, collapse = ", "))
  model_fitted <- mcglm(outcome, data = data,
                        family = family_outcome,
                        method = methods_fit,
                        p01 = p01, p10 = p10, pi_z = pi_z, Pi = Pi,
                        weights = case_weights)

  # ---- predictions and point estimate ---------------------------------
  y_nons_pred <- predict(model_fitted, type = "response",
                         method = mc_method)
  y_rand_pred <- predict(model_fitted, newdata = vars_rand,
                         type = "response", method = mc_method)
  residuals_y <- as.numeric(y_nons - y_nons_pred)

  svydesign_updated <- update(svydesign, .y_mc_hat = y_rand_pred)
  svydesign_mean    <- svymean(~.y_mc_hat, svydesign_updated)
  mu_hat <- as.numeric(svydesign_mean)

  # ---- variance --------------------------------------------------------
  if (se) {
    var_prob    <- as.numeric(attr(svydesign_mean, "var"))
    var_nonprob <- .var_nonprob_mc(model_fitted, mc_method, parts_rand,
                                   svydesign, family_outcome)
    var_total   <- var_prob + var_nonprob
    se_total    <- sqrt(var_total)
  } else {
    var_prob <- var_nonprob <- var_total <- se_total <- NA_real_
  }

  alpha <- 1 - 0.95
  ci <- mu_hat + c(-1, 1) * qnorm(1 - alpha / 2) * se_total

  # ---- assemble the nonprob-class object -------------------------------
  X_nons <- .build_xi(parts_nons$z_hat, parts_nons$x, parts_nons$K)
  X_rand <- .build_xi(parts_rand$z_hat, parts_rand$x, parts_nons$K)
  w_rand <- weights(svydesign)
  pop_size_fixed <- !is.null(pop_size)
  if (is.null(pop_size)) pop_size <- sum(w_rand)

  outcome_list <- list(
    model_fitted = model_fitted,
    coefficients = coef(model_fitted, method = mc_method),
    family       = family_outcome,
    mc_method    = mc_method,
    method       = "mcglm"
  )

  structure(
    list(
      call                = cl,
      data                = data,
      X                   = rbind(X_nons, X_rand),
      y                   = stats::setNames(list(y_nons), y_name),
      R                   = c(rep(1, n_nons), rep(0, n_rand)),
      ps_scores           = NULL,
      case_weights        = if (is.null(case_weights)) rep(1, n_nons)
                            else as.numeric(case_weights),
      ipw_weights         = NULL,
      control             = list(control_selection = NULL,
                                 control_outcome   = control_outcome,
                                 control_inference = control_inference),
      output              = data.frame(mean = mu_hat, SE = se_total,
                                       row.names = y_name),
      SE                  = data.frame(prob    = sqrt(var_prob),
                                       nonprob = sqrt(var_nonprob)),
      confidence_interval = data.frame(lower_bound = ci[1],
                                       upper_bound = ci[2],
                                       row.names = y_name),
      nonprob_size        = n_nons,
      prob_size           = n_rand,
      pop_size            = pop_size,
      pop_size_fixed      = pop_size_fixed,
      pop_totals          = NULL,
      pop_means           = NULL,
      outcome             = stats::setNames(list(outcome_list), y_name),
      selection           = NULL,
      boot_sample         = NULL,
      svydesign           = svydesign_updated,
      ys_rand_pred        = list(y_rand_pred),
      ys_nons_pred        = list(y_nons_pred),
      ys_resid            = list(residuals_y),
      estimator           = "mi",
      outcome_formula     = list(outcome),
      estimator_method    = paste0("mcglm (", family_outcome, ", ",
                                   mc_method, ")")
    ),
    class = c("nonprob_mc", "nonprob")
  )
}
