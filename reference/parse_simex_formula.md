# Parse a simex formula containing me() and mc() terms

Walks the formula parse tree, extracts measurement error /
misclassification descriptors, and produces a clean formula with
me()/mc() wrappers removed.

## Usage

``` r
parse_simex_formula(formula, data, env)
```

## Arguments

- formula:

  a formula potentially containing
  [`me()`](https://ojalab.github.io/mismeasured/reference/me.md) and
  [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) terms.

- data:

  the data frame.

- env:

  the calling environment (for evaluating sd / matrix arguments).

## Value

A list with components:

- clean_formula:

  A standard R formula with me/mc stripped.

- me_terms:

  List of me descriptors (each with `variable`, `sd` resolved to a
  numeric vector).

- mc_terms:

  List of mc descriptors (each with `variable`, `mc_matrix` resolved to
  a matrix).

- response_me:

  An me descriptor for the LHS, or NULL.

- error_type:

  Character: `"me"`, `"mc"`, or `"mixed"`.
