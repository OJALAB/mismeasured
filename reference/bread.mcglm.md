# Bread matrix for an `mcglm` fit

Returns the inverse of the (negated, averaged) Jacobian of the
estimating functions in
[`estfun.mcglm`](https://ojalab.github.io/mismeasured/reference/estfun.mcglm.md):
\\(\hat I + \hat M)^{-1}\\ for `method = "cs"` and \\\hat I^{-1}\\ for
`"naive"`, `"bca"`, `"bcm"`, each evaluated at that method's
coefficients.

## Usage

``` r
# S3 method for class 'mcglm'
bread(x, method = NULL, ...)
```

## Arguments

- x:

  An `mcglm` object.

- method:

  Estimation method; defaults to the last fitted one.

- ...:

  Unused.

## Value

A \\p \times p\\ matrix, such that `sandwich::sandwich(fit)` assembles
\\n^{-1} B \\\hat S\\ B^\top\\ with \\\hat S\\ from `estfun.mcglm`.

## See also

[`estfun.mcglm`](https://ojalab.github.io/mismeasured/reference/estfun.mcglm.md),
[`vcov.mcglm`](https://ojalab.github.io/mismeasured/reference/vcov.mcglm.md)
