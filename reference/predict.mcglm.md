# Linear-predictor or response predictions

Currently supports in-sample prediction only.

## Usage

``` r
# S3 method for class 'mcglm'
predict(
  object,
  newdata = NULL,
  type = c("link", "response"),
  method = NULL,
  ...
)
```

## Arguments

- object:

  An `mcglm` object.

- newdata:

  Optional new data. For formula fits, a data frame containing the
  [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md)
  variable (interpreted as the covariate at which to predict – pass the
  *true* category when available, e.g. when predicting onto a validation
  or probability sample) and the other covariates. For matrix-interface
  fits, a list with components `z_hat` (integer, coded `0, ..., K-1`)
  and `x` (covariate matrix as in the fit). `NULL` (default) predicts at
  the in-sample proxy design \\\hat\xi_i\\.

- type:

  One of `"link"`, `"response"`.

- method:

  Estimation method.

- ...:

  Unused.
