# Specify a misclassified variable

Use inside a
[`simex`](https://ojalab.github.io/mismeasured/reference/simex.md)
formula to mark a discrete covariate as subject to misclassification.

## Usage

``` r
mc(variable, matrix)
```

## Arguments

- variable:

  bare name of the factor variable in the data.

- matrix:

  a K x K misclassification matrix where
  `matrix[j, l] = P(observed = j | true = l)`. Columns must sum to 1.
  Can be a bare name of an object in the calling environment.

## Value

An object of class `"mc_term"` (used internally by
[`simex`](https://ojalab.github.io/mismeasured/reference/simex.md)).

## Category order and `matrix`

The rows and columns of `matrix` follow the *factor level order* of
`variable` (for integer-coded variables, the codes `0, ..., K-1`); the
first level is the baseline, category 0. This alignment cannot be
checked automatically: if you reorder the factor levels, you must
permute the rows and columns of `matrix` accordingly, otherwise the
correction is silently applied to the wrong categories.
[`summary()`](https://rdrr.io/r/base/summary.html) of an
[`mcglm`](https://ojalab.github.io/mismeasured/reference/mcglm.md) fit
prints the assumed category order for inspection.

## See also

[`me`](https://ojalab.github.io/mismeasured/reference/me.md),
[`simex`](https://ojalab.github.io/mismeasured/reference/simex.md)

## Examples

``` r
if (FALSE) { # \dontrun{
Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
simex(y ~ mc(z, Pi) + x, family = poisson(), data = df)
} # }
```
