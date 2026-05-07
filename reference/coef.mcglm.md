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
  `"bcm"`, `"cs"`, `"onestep"`). If `NULL`, returns a named list of all
  coefficient vectors.

- ...:

  Unused.
