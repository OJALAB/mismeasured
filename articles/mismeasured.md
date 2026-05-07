# Getting started with mismeasured

## What this package does

**mismeasured** corrects regression estimates for measurement error and
misclassification of covariates (and response). It bundles two
complementary toolsets:

| Toolset | When to use | Vignette |
|----|----|----|
| [`simex()`](https://ojalab.github.io/mismeasured/reference/simex.md) | Continuous measurement error, Berkson error, single or multiple misclassified factors, response misclassification | [`vignette("simex", "mismeasured")`](https://ojalab.github.io/mismeasured/articles/simex.md) |
| [`mcglm()`](https://ojalab.github.io/mismeasured/reference/mcglm.md) / [`mclm()`](https://ojalab.github.io/mismeasured/reference/mclm.md) | Misclassified discrete covariate(s) with known (or auxiliary-estimated) $`\Pi`$, where you want analytical bias-corrected GLM/LM inference | [`vignette("mcglm", "mismeasured")`](https://ojalab.github.io/mismeasured/articles/mcglm.md) |

Variables subject to error are marked directly in the formula:
`me(x, sd)` for continuous measurement error, `mc(z, Pi)` for discrete
misclassification.

### Supported response families

[`gaussian()`](https://rdrr.io/r/stats/family.html),
[`poisson()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html) (canonical links).
[`mcglm()`](https://ojalab.github.io/mismeasured/reference/mcglm.md)
also supports `family = "multinomial"` for the naive and one-step
estimators.

## A SIMEX example in 30 seconds

``` r

n <- 1500
x_true <- rnorm(n)
y      <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
x_obs  <- x_true + rnorm(n, sd = 0.5)
df     <- data.frame(y = y, x = x_obs)

fit_simex <- simex(y ~ me(x, 0.5), data = df, B = 200, seed = 42)

c(naive_lm = coef(lm(y ~ x, data = df))["x"],
  simex    = coef(fit_simex)["x"],
  truth    = 2)
#> naive_lm.x    simex.x      truth 
#>   1.578317   1.921585   2.000000
```

The naive OLS slope is attenuated; SIMEX recovers a value close to the
true slope of 2. See
[`vignette("simex", "mismeasured")`](https://ojalab.github.io/mismeasured/articles/simex.md)
for the full algorithm, MC-SIMEX, response misclassification, Berkson
error, and weights.

## A bias-corrected GLM example in 30 seconds

``` r

n  <- 3000
Pi <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
z      <- rbinom(n, 1, 0.4)
z_hat  <- z
z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.10)
x1 <- rnorm(n)
y  <- rpois(n, exp(0.8 * z - 0.5 + 0.7 * x1))
df <- data.frame(y = y, z = factor(z_hat), x1 = x1)

fit_mcglm <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
                   method = c("naive", "cs", "cs_akn"))
coef(fit_mcglm, method = "naive")
#>      gamma     alpha0     alpha1 
#>  0.6201552 -0.4452069  0.7125682
coef(fit_mcglm, method = "cs")
#>      gamma     alpha0     alpha1 
#>  0.7792699 -0.5291241  0.7141844
coef(fit_mcglm, method = "cs_akn")
#>      gamma     alpha0     alpha1 
#>  0.7789926 -0.5294969  0.7161404
```

The naive $`\hat\gamma`$ is attenuated below the truth (0.8); both
corrected scores recover it.
[`mclm()`](https://ojalab.github.io/mismeasured/reference/mclm.md) is a
Gaussian-fixed wrapper of
[`mcglm()`](https://ojalab.github.io/mismeasured/reference/mcglm.md) for
linear models. See
[`vignette("mcglm", "mismeasured")`](https://ojalab.github.io/mismeasured/articles/mcglm.md)
for the full theory, the BCA/BCM/one-step alternatives, and the
multicategory case.

## Utility functions

The package also exports two helpers for working with misclassification
matrices, used by both toolsets:

- [`build.mc.matrix()`](https://ojalab.github.io/mismeasured/reference/build.mc.matrix.md)
  finds the nearest valid (non-negative-power) matrix given an empirical
  estimate of $`\Pi`$.
- [`check.mc.matrix()`](https://ojalab.github.io/mismeasured/reference/check.mc.matrix.md)
  validates a list of candidate $`\Pi`$s.

``` r

Pi_empirical <- matrix(c(0.85, 0.10, 0.05,
                         0.10, 0.80, 0.10,
                         0.05, 0.10, 0.85), 3, 3, byrow = TRUE)
check.mc.matrix(list(Pi_empirical))
#> [1] TRUE
Pi_valid <- build.mc.matrix(Pi_empirical)
check.mc.matrix(list(Pi_valid))
#> [1] TRUE
```

## Getting help

- [`?simex`](https://ojalab.github.io/mismeasured/reference/simex.md),
  [`?mcglm`](https://ojalab.github.io/mismeasured/reference/mcglm.md),
  [`?mclm`](https://ojalab.github.io/mismeasured/reference/mclm.md),
  [`?me`](https://ojalab.github.io/mismeasured/reference/me.md),
  [`?mc`](https://ojalab.github.io/mismeasured/reference/mc.md) —
  function references.
- [`vignette("simex", "mismeasured")`](https://ojalab.github.io/mismeasured/articles/simex.md)
  — SIMEX/MC-SIMEX in depth.
- [`vignette("mcglm", "mismeasured")`](https://ojalab.github.io/mismeasured/articles/mcglm.md)
  — bias-corrected GLM theory and methods.
- `https://ojalab.github.io/mismeasured/` — pkgdown site with the full
  reference index.
- Issues: `https://github.com/OJALAB/mismeasured/issues`.
