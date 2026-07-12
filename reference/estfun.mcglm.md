# Extract empirical estimating functions from an `mcglm` fit

Returns the \\n \times p\\ matrix of per-observation estimating
functions evaluated at the coefficients of the chosen `method`. For
`method = "cs"` these are the corrected-score contributions
\\\phi_i(\hat\psi) = \hat\xi_i\\y_i - \mu(\tilde\eta_i)\\ -
m_i(\hat\psi)\\; for `"naive"`, `"bca"` and `"bcm"` they are the
proxy-score contributions \\\hat\xi_i\\y_i - \mu(\tilde\eta_i)\\\\
(whose sandwich is the asymptotic variance of all three estimators in
the drifting regime).

## Usage

``` r
# S3 method for class 'mcglm'
estfun(x, method = NULL, ...)
```

## Arguments

- x:

  An `mcglm` object.

- method:

  Estimation method; defaults to the last fitted one. One of `"naive"`,
  `"bca"`, `"bcm"`, `"cs"`.

- ...:

  Unused.

## Value

An \\n \times p\\ numeric matrix. With frequency weights the rows are
multiplied by the weights (the sandwich convention for `glm`), in which
case
[`sandwich::sandwich()`](https://sandwich.R-Forge.R-project.org/reference/sandwich.html)
treats them as probability weights and differs from
[`vcov()`](https://rdrr.io/r/stats/vcov.html), which treats them as
frequencies.

## Details

Together with
[`bread.mcglm`](https://ojalab.github.io/mismeasured/reference/bread.mcglm.md)
this makes `mcglm` fits work with the sandwich ecosystem (`meat()`,
`lmtest::coeftest(vcov. = )`, bootstrap tools, ...).

## Assembling the sandwich

Prefer `vcov(fit, method = )`: it already *is* the correct sandwich
\\J^{-1} S J^{-\top}/n\\. Do not pass `method` to
[`sandwich::sandwich()`](https://sandwich.R-Forge.R-project.org/reference/sandwich.html)
– it forwards extra arguments to the meat but *not* to the bread,
silently mixing estimators – and note that
[`sandwich::sandwich()`](https://sandwich.R-Forge.R-project.org/reference/sandwich.html)
assembles `bread %*% meat %*% bread` without a transpose, which is only
valid for the symmetric bread of `"naive"`/`"bca"`/`"bcm"`; the
corrected-score bread \\(\hat I + \hat M)^{-1}\\ is not symmetric. To
assemble by hand: `B %*% meat(fit, method = m) %*% t(B) / n` with
`B = bread(fit, method = m)`, which reproduces `vcov(fit, method = m)`
exactly for unweighted fits.

## See also

[`bread.mcglm`](https://ojalab.github.io/mismeasured/reference/bread.mcglm.md),
[`vcov.mcglm`](https://ojalab.github.io/mismeasured/reference/vcov.mcglm.md)
