# Extract coefficients from an `mcglm` fit

Extract coefficients from an `mcglm` fit

## Usage

``` r
# S3 method for class 'mcglm'
coef(object, method = NULL, ...)
```

## Arguments

- object:

  An `mcglm` object.

- method:

  Optional character: which estimator to extract (`"naive"`, `"bca"`,
  `"bcm"`, `"cs"`, `"onestep"`). If `NULL`, returns the coefficients of
  the default method (the last one fit), consistently with
  [`vcov.mcglm`](https://ojalab.github.io/mismeasured/reference/vcov.mcglm.md)
  and
  [`predict.mcglm`](https://ojalab.github.io/mismeasured/reference/predict.mcglm.md)
  – and as required by downstream tools such as
  [`lmtest::coeftest`](https://rdrr.io/pkg/lmtest/man/coeftest.html) and
  marginaleffects. Use `method = "all"` for the named list of all fitted
  coefficient vectors (the pre-0.6.0 default).

- ...:

  Unused.
