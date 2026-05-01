# =========================================================================
# Misclassification utilities
# =========================================================================

#' Generate misclassified data
#'
#' Takes a data frame of factor variables and produces misclassified versions
#' using the misclassification matrix raised to power \code{k} via
#' eigendecomposition.
#'
#' @param data.org data frame containing factor variables.
#' @param mc.matrix a named list of misclassification matrices. Names must
#'   correspond to variable names in \code{data.org}. Column names must match
#'   factor levels.
#' @param k exponent for the misclassification matrix (default: 1).
#'
#' @return A data frame with the misclassified variables.
#'
#' @examples
#' x1 <- factor(rbinom(100, 1, 0.5))
#' p1 <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2)
#' colnames(p1) <- levels(x1)
#' x <- data.frame(x1 = x1)
#' x.mc <- misclass(x, list(x1 = p1), k = 1)
#'
#' @export
misclass <- function(data.org, mc.matrix, k = 1) {
  if (!is.list(mc.matrix))
    stop("mc.matrix must be a list", call. = FALSE)
  if (!is.data.frame(data.org))
    stop("data.org must be a data.frame", call. = FALSE)
  if (!all(names(mc.matrix) %in% colnames(data.org)))
    stop("Names of mc.matrix and colnames of data.org do not match",
         call. = FALSE)
  if (k < 0)
    stop("k must be non-negative")

  data.mc <- data.org
  factors <- lapply(data.org, levels)
  ev <- lapply(mc.matrix, eigen)

  for (j in names(mc.matrix)) {
    evalue <- ev[[j]]$values
    evectors <- ev[[j]]$vectors
    d <- diag(evalue)
    mc <- zapsmall(evectors %*% d^k %*% solve(evectors))
    dimnames(mc) <- dimnames(mc.matrix[[j]])
    for (i in factors[[j]]) {
      data.mc[[j]][data.org[[j]] == i] <- sample(
        x = factors[[j]],
        size = sum(data.org[[j]] == i),
        prob = mc[, i],
        replace = TRUE
      )
    }
  }
  data.mc
}


#' Build a valid misclassification matrix
#'
#' Estimates the nearest misclassification matrix that can be raised to
#' fractional powers (i.e., has a valid matrix logarithm/generator).
#'
#' @param mc.matrix an empirical misclassification matrix.
#' @param method method for estimating the generator: \code{"series"} (default),
#'   \code{"log"}, or \code{"jlt"} (Jarrow-Lando-Turnbull).
#' @param tuning small positive constant for numerical stability (default:
#'   \code{sqrt(.Machine$double.eps)}).
#' @param diag.cor logical; if \code{TRUE}, corrections are subtracted from the
#'   diagonal only. If \code{FALSE} (default), distributed proportionally.
#' @param tol convergence tolerance for the series method.
#' @param max.iter maximum iterations for the series method.
#'
#' @return A valid misclassification matrix.
#'
#' @references
#' Israel, R.B., Rosenthal, J.S., Wei, J.Z. (2001). Finding generators for
#' Markov Chains via empirical transition matrices, with applications to credit
#' ratings. \emph{Mathematical Finance}, 11, 245--265.
#'
#' @examples
#' Pi <- matrix(c(0.989, 0.01, 0.001, 0.17, 0.829, 0.001, 0.001, 0.18, 0.819),
#'              nrow = 3, byrow = FALSE)
#' check.mc.matrix(list(Pi))
#' Pi_fixed <- build.mc.matrix(Pi)
#' check.mc.matrix(list(Pi_fixed))
#'
#' @export
build.mc.matrix <- function(mc.matrix, method = "series",
                            tuning = sqrt(.Machine$double.eps),
                            diag.cor = FALSE,
                            tol = .Machine$double.eps,
                            max.iter = 100) {
  nam <- dimnames(mc.matrix)

  if (method == "jlt") {
    if (any(diag(mc.matrix) == 0))
      stop("0 on the diagonal not allowed for method 'jlt'")
    generator <- t(t(mc.matrix) * log(diag(mc.matrix)) /
                      (diag(mc.matrix) - 1))
    diag(generator) <- log(diag(mc.matrix))
  }

  if (method == "series") {
    iter <- 1
    eps <- 1
    ident <- diag(nrow(mc.matrix))
    generator <- matrix(0, nrow = nrow(mc.matrix), ncol = ncol(mc.matrix))
    while (sum(abs(eps)) > tol && iter < max.iter) {
      eps <- ident
      for (i in seq_len(iter)) eps <- eps %*% (mc.matrix - ident)
      eps <- (-1)^(iter + 1) * eps / iter
      generator <- generator + eps
      iter <- iter + 1
    }
    if (iter == max.iter)
      stop("Series did not converge, try method = 'jlt'")
  }

  if (method == "log") {
    ev <- eigen(mc.matrix)
    generator <- ev$vectors %*% diag(log(ev$values)) %*% solve(ev$vectors)
  }

  # Correct negative off-diagonal entries in the generator
  if (diag.cor) {
    for (i in seq_len(nrow(mc.matrix))) {
      if (any(generator[-i, i] < 0)) {
        idx <- generator[, i] < 0
        corrector <- rep(0, nrow(generator))
        corrector[idx] <- abs(generator[idx, i]) + tuning
        idx[i] <- FALSE
        corrector[i] <- 0
        generator[, i] <- generator[, i] + corrector
        generator[i, i] <- generator[i, i] - sum(corrector)
      }
    }
  } else {
    for (i in seq_len(nrow(mc.matrix))) {
      corrector <- rep(0, nrow(generator))
      g_val <- abs(generator[i, i]) + sum(pmax(generator[-i, i], 0))
      b_val <- sum(pmax(-generator[-i, i] + tuning, 0))
      idx <- generator[, i] < 0
      corrector[idx] <- tuning
      corrector[i] <- generator[i, i]
      idx2 <- g_val > 0
      idx[i] <- FALSE
      if (g_val > 0) {
        mask <- !idx & (seq_len(nrow(generator)) != i)
        corrector[mask] <- (generator[mask, i] -
                              (b_val * abs(generator[mask, i]) / g_val))
      }
      generator[, i] <- corrector
    }
  }

  # Reconstruct the matrix via expm(generator)
  ev <- eigen(generator)
  mc.matrix <- zapsmall(ev$vectors %*% diag(exp(ev$values)) %*%
                          solve(ev$vectors))
  dimnames(mc.matrix) <- nam
  mc.matrix
}


#' Check if a misclassification matrix is valid for fractional powers
#'
#' Tests whether the misclassification matrix can be raised to powers less
#' than 1 without producing negative probabilities.
#'
#' @param mc.matrix a list of misclassification matrices.
#' @param tol tolerance for checking non-negativity (default:
#'   \code{.Machine$double.eps}).
#'
#' @return A logical vector indicating validity for each matrix.
#'
#' @examples
#' P1 <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2)
#' P2 <- matrix(c(0.4, 0.6, 0.6, 0.4), nrow = 2)
#' check.mc.matrix(list(P1, P2))  # TRUE FALSE
#'
#' @export
check.mc.matrix <- function(mc.matrix, tol = .Machine$double.eps) {
  result <- logical(length(mc.matrix))
  for (i in seq_along(mc.matrix)) {
    m <- mc.matrix[[i]]
    if (all(dim(m) == c(2, 2))) {
      result[i] <- (m[1, 1] + m[2, 2] > 1)
    } else {
      ev <- eigen(m)
      mc_log <- zapsmall(ev$vectors %*% diag(log(ev$values)) %*%
                           solve(ev$vectors))
      diag(mc_log) <- 1
      result[i] <- all(mc_log > -tol)
    }
  }
  result
}


#' Construct a block diagonal matrix
#'
#' @param d a list of matrices/vectors, or a single matrix/vector.
#' @param n number of repetitions (used when \code{d} is not a list).
#'
#' @return A block diagonal matrix.
#'
#' @examples
#' a <- matrix(1, 2, 2)
#' b <- matrix(2, 2, 3)
#' diag.block(list(a, b))
#'
#' @export
diag.block <- function(d, n) {
  if (is.list(d)) {
    d.row <- sapply(d, NROW)
    d.col <- sapply(d, NCOL)
    result <- matrix(0, nrow = sum(d.row), ncol = sum(d.col))
    d.row <- c(0, cumsum(d.row))
    d.col <- c(0, cumsum(d.col))
    for (i in seq_along(d)) {
      result[(d.row[i] + 1):d.row[i + 1],
             (d.col[i] + 1):d.col[i + 1]] <- as.matrix(d[[i]])
    }
  } else {
    result <- kronecker(diag(1, n), d)
  }
  result
}
