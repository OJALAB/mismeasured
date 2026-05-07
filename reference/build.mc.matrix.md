# Build a valid misclassification matrix

Estimates the nearest misclassification matrix that can be raised to
fractional powers (i.e., has a valid matrix logarithm/generator).

## Usage

``` r
build.mc.matrix(
  mc.matrix,
  method = "series",
  tuning = sqrt(.Machine$double.eps),
  diag.cor = FALSE,
  tol = .Machine$double.eps,
  max.iter = 100
)
```

## Arguments

- mc.matrix:

  an empirical misclassification matrix.

- method:

  method for estimating the generator: `"series"` (default), `"log"`, or
  `"jlt"` (Jarrow-Lando-Turnbull).

- tuning:

  small positive constant for numerical stability (default:
  `sqrt(.Machine$double.eps)`).

- diag.cor:

  logical; if `TRUE`, corrections are subtracted from the diagonal only.
  If `FALSE` (default), distributed proportionally.

- tol:

  convergence tolerance for the series method.

- max.iter:

  maximum iterations for the series method.

## Value

A valid misclassification matrix.

## References

Israel, R.B., Rosenthal, J.S., Wei, J.Z. (2001). Finding generators for
Markov Chains via empirical transition matrices, with applications to
credit ratings. *Mathematical Finance*, 11, 245–265.

## Examples

``` r
Pi <- matrix(c(0.989, 0.01, 0.001, 0.17, 0.829, 0.001, 0.001, 0.18, 0.819),
             nrow = 3, byrow = FALSE)
check.mc.matrix(list(Pi))
#> [1] FALSE
Pi_fixed <- build.mc.matrix(Pi)
check.mc.matrix(list(Pi_fixed))
#> [1] FALSE
```
