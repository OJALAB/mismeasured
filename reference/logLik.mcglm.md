# Log-likelihood of an `mcglm` fit

Returns the naive-GLM log-likelihood for `method = "naive"` (the only
estimator with a closed-form likelihood for proxy data) or the
integrated mixture log-likelihood for `method = "onestep"`. Other
methods do not correspond to a likelihood; an error is thrown.

## Usage

``` r
# S3 method for class 'mcglm'
logLik(object, method = NULL, ...)
```

## Arguments

- object:

  An `mcglm` object.

- method:

  Estimation method.

- ...:

  Unused.
