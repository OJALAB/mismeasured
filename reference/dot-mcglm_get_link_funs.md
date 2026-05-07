# Get inverse-link, its derivative, and second derivative for a GLM family

Get inverse-link, its derivative, and second derivative for a GLM family

## Usage

``` r
.mcglm_get_link_funs(family)
```

## Arguments

- family:

  A [`family`](https://rdrr.io/r/stats/family.html) object or character
  string ("poisson", "binomial", "gaussian").

## Value

A list with components `mu` (inverse link), `mu_dot` (first derivative),
`mu_ddot` (second derivative), and `family`.
