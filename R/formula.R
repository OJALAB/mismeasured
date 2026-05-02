# =========================================================================
# Formula parsing infrastructure for SIMEX
# =========================================================================
# Provides me() and mc() formula terms inspired by brms::me().
# parse_simex_formula() walks the formula tree, extracts descriptors,
# and returns a clean glm-compatible formula.

#' Specify a variable measured with error
#'
#' Use inside a \code{\link{simex}} formula to mark a covariate (or response)
#' as measured with additive Gaussian error.
#'
#' @param variable bare name of the variable in the data.
#' @param sd measurement error standard deviation. Can be a scalar
#'   (homoscedastic), a numeric vector of length n, or a bare column name
#'   in the data (heteroscedastic).
#' @param type type of measurement error: \code{"classical"} (default) for
#'   the standard additive error model \eqn{W = X + U}, or \code{"berkson"}
#'   for the Berkson error model \eqn{X = W + U} where the true value is
#'   a noisy version of the observed.
#' @param mean mean of the measurement error distribution (default 0).
#'   Use for systematic (non-zero mean) measurement error.
#'
#' @return An object of class \code{"me_term"} (used internally by
#'   \code{\link{simex}}; not intended to be called directly).
#'
#' @examples
#' \dontrun{
#' simex(y ~ me(x, 0.5) + w, data = df)
#' simex(y ~ me(x, sd_x) + w, data = df)            # heteroscedastic
#' simex(y ~ me(x, 0.5, type = "berkson") + w, data = df)  # Berkson error
#' simex(y ~ me(x, 0.5, mean = 0.1) + w, data = df)       # non-zero mean
#' }
#'
#' @seealso \code{\link{mc}}, \code{\link{simex}}
#' @export
me <- function(variable, sd, type = "classical", mean = 0) {
  type <- match.arg(type, c("classical", "berkson"))
  structure(
    list(
      variable = deparse(substitute(variable)),
      sd_expr  = substitute(sd),
      type     = type,
      mean     = mean
    ),
    class = "me_term"
  )
}


#' Specify a misclassified variable
#'
#' Use inside a \code{\link{simex}} formula to mark a discrete covariate
#' as subject to misclassification.
#'
#' @param variable bare name of the factor variable in the data.
#' @param matrix a K x K misclassification matrix where
#'   \code{matrix[j, l] = P(observed = j | true = l)}. Columns must sum to 1.
#'   Can be a bare name of an object in the calling environment.
#'
#' @return An object of class \code{"mc_term"} (used internally by
#'   \code{\link{simex}}).
#'
#' @examples
#' \dontrun{
#' Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
#' simex(y ~ mc(z, Pi) + x, family = poisson(), data = df)
#' }
#'
#' @seealso \code{\link{me}}, \code{\link{simex}}
#' @export
mc <- function(variable, matrix) {
  structure(
    list(
      variable = deparse(substitute(variable)),
      mat_expr = substitute(matrix)
    ),
    class = "mc_term"
  )
}


# =========================================================================
# Formula parser
# =========================================================================

#' Parse a simex formula containing me() and mc() terms
#'
#' Walks the formula parse tree, extracts measurement error / misclassification
#' descriptors, and produces a clean formula with me()/mc() wrappers removed.
#'
#' @param formula a formula potentially containing \code{me()} and \code{mc()}
#'   terms.
#' @param data the data frame.
#' @param env the calling environment (for evaluating sd / matrix arguments).
#'
#' @return A list with components:
#'   \describe{
#'     \item{clean_formula}{A standard R formula with me/mc stripped.}
#'     \item{me_terms}{List of me descriptors (each with \code{variable},
#'       \code{sd} resolved to a numeric vector).}
#'     \item{mc_terms}{List of mc descriptors (each with \code{variable},
#'       \code{mc_matrix} resolved to a matrix).}
#'     \item{response_me}{An me descriptor for the LHS, or NULL.}
#'     \item{error_type}{Character: \code{"me"}, \code{"mc"}, or \code{"mixed"}.}
#'   }
#'
#' @keywords internal
parse_simex_formula <- function(formula, data, env) {
  # Separate LHS and RHS
  lhs <- formula[[2]]
  rhs <- formula[[3]]

  # --- Walk LHS for me() or mc() ---
  response_me <- NULL
  response_mc <- NULL
  if (is.call(lhs) && identical(lhs[[1]], as.name("me"))) {
    var_name <- deparse(lhs[[2]])
    sd_val <- .resolve_sd(lhs[[3]], var_name, data, env)
    response_me <- list(variable = var_name, sd = sd_val)
    lhs <- lhs[[2]]  # strip me(), keep bare variable name
    stop("Response measurement error (me() on the left-hand side) is not yet supported.",
         call. = FALSE)
  }
  if (is.call(lhs) && identical(lhs[[1]], as.name("mc"))) {
    var_name <- deparse(lhs[[2]])
    mat_val <- .resolve_matrix(lhs[[3]], data, env)
    # Validate response mc matrix
    if (!is.factor(data[[var_name]])) {
      warning("mc() response variable '", var_name,
              "' is not a factor; coercing.", call. = FALSE)
      data[[var_name]] <- factor(data[[var_name]])
    }
    K_resp <- nrow(mat_val)
    n_levels_resp <- length(levels(data[[var_name]]))
    if (K_resp != n_levels_resp)
      stop("mc() response matrix has ", K_resp, " rows but variable '",
           var_name, "' has ", n_levels_resp, " levels.", call. = FALSE)
    cs <- colSums(mat_val)
    if (any(abs(cs - 1) > 1e-6))
      stop("Columns of mc() response matrix must sum to 1 (got: ",
           paste(round(cs, 4), collapse = ", "), ").", call. = FALSE)
    response_mc <- list(variable = var_name, mc_matrix = mat_val)
    lhs <- lhs[[2]]  # strip mc(), keep bare variable name
  }

  # --- Walk RHS for me() and mc() ---
  descriptors <- list()
  rhs_clean <- .walk_and_collect(rhs, data, env, descriptors = descriptors,
                                  depth = 0, inside_interaction = FALSE)

  # Collect results from the walk
  me_terms <- Filter(function(d) d$type == "me", .env_descriptors$items)
  mc_terms <- Filter(function(d) d$type == "mc", .env_descriptors$items)

  # --- Validation ---

  for (mt in me_terms) {
    if (is.factor(data[[mt$variable]]))
      stop("me() variable '", mt$variable,
           "' is a factor. Use mc() for discrete misclassification.", call. = FALSE)
  }

  for (mt in mc_terms) {
    if (!is.factor(data[[mt$variable]])) {
      warning("mc() variable '", mt$variable,
              "' is not a factor; coercing.", call. = FALSE)
      data[[mt$variable]] <- factor(data[[mt$variable]])
    }
    K <- nrow(mt$mc_matrix)
    n_levels <- length(levels(data[[mt$variable]]))
    if (K != n_levels)
      stop("mc() matrix has ", K, " rows but variable '", mt$variable,
           "' has ", n_levels, " levels.", call. = FALSE)
    cs <- colSums(mt$mc_matrix)
    if (any(abs(cs - 1) > 1e-6))
      stop("Columns of mc() matrix must sum to 1 (got: ",
           paste(round(cs, 4), collapse = ", "), ").", call. = FALSE)
  }

  # --- Determine error type ---
  has_me <- length(me_terms) > 0
  has_mc <- length(mc_terms) > 0 || !is.null(response_mc)
  if (has_me && has_mc) {
    stop("Mixed me() and mc() terms in the same formula are not yet supported.",
         call. = FALSE)
  }
  if (!has_me && !has_mc) {
    stop("Formula contains no me() or mc() terms. ",
         "Use glm() for models without measurement error.", call. = FALSE)
  }
  error_type <- if (has_me) "me" else "mc"

  # Rebuild clean formula
  clean_formula <- as.formula(call("~", lhs, rhs_clean), env = env)

  list(
    clean_formula = clean_formula,
    me_terms      = me_terms,
    mc_terms      = mc_terms,
    response_me   = response_me,
    response_mc   = response_mc,
    error_type    = error_type
  )
}


# =========================================================================
# Formula tree walker
# =========================================================================

# We use a simple environment to accumulate descriptors during the walk,
# since R's recursive functions don't have mutable state otherwise.

#' Walk a formula node, collect me/mc descriptors, return cleaned node
#' @keywords internal
.walk_and_collect <- function(node, data, env, descriptors, depth,
                               inside_interaction) {
  # Initialize accumulator on first call
  if (depth == 0) {
    .env_descriptors$items <- list()
  }

  if (is.name(node) || is.numeric(node) || is.character(node)) {
    return(node)
  }

  if (!is.call(node)) return(node)

  fn_name <- deparse(node[[1]])

  # --- Detect me() ---
  if (fn_name == "me") {
    if (inside_interaction)
      stop("me() terms cannot appear inside interactions (: or *).",
           call. = FALSE)
    if (length(node) < 3)
      stop("me() requires two arguments: me(variable, sd).", call. = FALSE)

    var_name <- deparse(node[[2]])
    sd_val <- .resolve_sd(node[[3]], var_name, data, env)

    # Extract optional type and mean arguments
    me_call <- match.call(me, node)
    me_type <- if (!is.null(me_call$type))
      eval(me_call$type, envir = data, enclos = env) else "classical"
    me_mean <- if (!is.null(me_call$mean))
      eval(me_call$mean, envir = data, enclos = env) else 0

    .env_descriptors$items <- c(.env_descriptors$items, list(
      list(type = "me", variable = var_name, sd = sd_val,
           error_type = me_type, mean = me_mean)
    ))

    return(node[[2]])  # replace me(x, sd) with x
  }

  # --- Detect mc() ---
  if (fn_name == "mc") {
    if (inside_interaction)
      stop("mc() terms cannot appear inside interactions (: or *).",
           call. = FALSE)
    if (length(node) < 3)
      stop("mc() requires two arguments: mc(variable, matrix).", call. = FALSE)

    var_name <- deparse(node[[2]])
    mat_val <- .resolve_matrix(node[[3]], data, env)

    .env_descriptors$items <- c(.env_descriptors$items, list(
      list(type = "mc", variable = var_name, mc_matrix = mat_val)
    ))

    return(node[[2]])  # replace mc(z, Pi) with z
  }

  # --- Detect interaction operators ---
  is_interact <- fn_name %in% c(":", "*")

  # --- Detect wrapping functions like I() ---
  if (fn_name == "I") {
    # Check if me/mc is inside I() — not allowed
    inner <- node[[2]]
    if (is.call(inner) && deparse(inner[[1]]) %in% c("me", "mc"))
      stop("me()/mc() terms cannot appear inside I().", call. = FALSE)
  }

  # --- Recurse into sub-nodes ---
  for (i in seq_along(node)[-1]) {
    node[[i]] <- .walk_and_collect(node[[i]], data, env, descriptors,
                                    depth + 1,
                                    inside_interaction || is_interact)
  }

  node
}

# Accumulator environment (package-level, reset on each parse)
.env_descriptors <- new.env(parent = emptyenv())
.env_descriptors$items <- list()


# =========================================================================
# Argument resolution helpers
# =========================================================================

#' Resolve the sd argument of me()
#'
#' Evaluates in data first (column name), then in calling env.
#' @keywords internal
.resolve_sd <- function(sd_expr, var_name, data, env) {
  # Try as column name in data
  sd_name <- tryCatch(deparse(sd_expr), error = function(e) NULL)
  if (!is.null(sd_name) && sd_name %in% names(data)) {
    sd_val <- data[[sd_name]]
    if (!is.numeric(sd_val))
      stop("me() sd column '", sd_name, "' must be numeric.", call. = FALSE)
    return(as.numeric(sd_val))
  }

  # Evaluate in calling environment
  sd_val <- tryCatch(
    eval(sd_expr, envir = data, enclos = env),
    error = function(e) {
      stop("Cannot resolve me() sd argument: ", deparse(sd_expr),
           call. = FALSE)
    }
  )

  if (!is.numeric(sd_val))
    stop("me() sd must be numeric.", call. = FALSE)
  if (any(sd_val < 0))
    stop("me() sd must be non-negative.", call. = FALSE)

  as.numeric(sd_val)
}


#' Resolve the matrix argument of mc()
#' @keywords internal
.resolve_matrix <- function(mat_expr, data, env) {
  mat_val <- tryCatch(
    eval(mat_expr, envir = data, enclos = env),
    error = function(e) {
      stop("Cannot resolve mc() matrix argument: ", deparse(mat_expr),
           call. = FALSE)
    }
  )

  if (!is.matrix(mat_val))
    mat_val <- as.matrix(mat_val)

  if (nrow(mat_val) != ncol(mat_val))
    stop("mc() matrix must be square (K x K).", call. = FALSE)

  mat_val
}
