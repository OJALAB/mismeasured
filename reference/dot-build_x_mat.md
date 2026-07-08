# Build x_mat from a clean formula, preserving I(), interactions, offsets, etc.

Replaces the fragile all.vars() + reformulate() approach. Uses
model.matrix() on the clean formula and removes columns belonging to mc
variables.

## Usage

``` r
.build_x_mat(clean_formula, mc_variables, data, xlev = NULL)
```

## Arguments

- clean_formula:

  the formula with mc() wrappers stripped (bare variable names remain)

- mc_variables:

  character vector of mc variable names to exclude

- data:

  the data frame

- xlev:

  a named list of factor levels to enforce (as in
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)); used by
  `predict.simex` so newdata is coded with the training levels

## Value

a list with x_mat (the non-mc part of the model matrix) and x_terms (the
terms object, for use in predict)
