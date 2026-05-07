# Check if a misclassification matrix is valid for fractional powers

Tests whether the misclassification matrix can be raised to powers less
than 1 without producing negative probabilities.

## Usage

``` r
check.mc.matrix(mc.matrix, tol = .Machine$double.eps)
```

## Arguments

- mc.matrix:

  a list of misclassification matrices.

- tol:

  tolerance for checking non-negativity (default:
  `.Machine$double.eps`).

## Value

A logical vector indicating validity for each matrix.

## Examples

``` r
P1 <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2)
P2 <- matrix(c(0.4, 0.6, 0.6, 0.4), nrow = 2)
check.mc.matrix(list(P1, P2))  # TRUE FALSE
#> [1]  TRUE FALSE
```
