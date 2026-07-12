# Parse an `mc()` model formula

Public interface to the formula parsing used by
[`mcglm`](https://ojalab.github.io/mismeasured/reference/mcglm.md):
splits `y ~ mc(z, Pi) + x1 + ...` into the response, the 0-based proxy
covariate, the model matrix of the remaining terms, and the
misclassification matrix. Intended for packages building on mismeasured
(e.g. survey-integration wrappers) that must intercept the
[`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) term
before standard formula machinery (such as
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)) sees it.

## Usage

``` r
mc_parse_formula(formula, data, env = parent.frame())
```

## Arguments

- formula:

  A formula with exactly one
  [`mc`](https://ojalab.github.io/mismeasured/reference/mc.md) term on
  the right-hand side, e.g. `y ~ mc(z, Pi) + x1`.

- data:

  A data frame containing the variables.

- env:

  Environment in which to evaluate the
  [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) matrix
  argument (default: the caller's environment).

## Value

A list with components `y` (numeric response), `z_hat` (integer proxy,
coded `0, ..., K-1`), `x` (model matrix of the remaining terms,
including an intercept), `Pi` (the \\K \times K\\ misclassification
matrix, or `NULL` if not supplied inside
[`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md)), and `K`
(number of categories).

## See also

[`mc`](https://ojalab.github.io/mismeasured/reference/mc.md),
[`mcglm`](https://ojalab.github.io/mismeasured/reference/mcglm.md)

## Examples

``` r
Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
df <- data.frame(y = rpois(50, 1), z = rbinom(50, 1, 0.4),
                 x1 = rnorm(50))
parts <- mc_parse_formula(y ~ mc(z, Pi) + x1, df)
str(parts)
#> List of 5
#>  $ y    : num [1:50] 0 2 1 0 0 1 1 0 1 2 ...
#>  $ z_hat: int [1:50] 0 0 0 0 0 1 0 1 0 0 ...
#>  $ x    : num [1:50, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:50] "1" "2" "3" "4" ...
#>   .. ..$ : chr [1:2] "(Intercept)" "x1"
#>   ..- attr(*, "assign")= int [1:2] 0 1
#>  $ Pi   : num [1:2, 1:2] 0.9 0.1 0.15 0.85
#>  $ K    : int 2
```
