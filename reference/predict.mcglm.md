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

  Currently unused. Predictions are computed at the in-sample proxy
  design \\\hat\xi_i\\.

- type:

  One of `"link"`, `"response"`.

- method:

  Estimation method.

- ...:

  Unused.
