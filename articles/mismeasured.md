# Correcting measurement error and misclassification with mismeasured

## Introduction

The **mismeasured** package provides Simulation-Extrapolation (SIMEX)
and Misclassification SIMEX (MC-SIMEX) corrections for generalized
linear models where covariates are measured with error or subject to
misclassification.

The package uses a formula interface inspired by **brms**: error-prone
variables are marked directly in the formula with
[`me()`](https://ojalab.github.io/mismeasured/reference/me.md)
(continuous measurement error) or
[`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) (discrete
misclassification). The simulation step runs in C++ via Rcpp and
RcppEigen for high performance.

### Supported error types

| Type | Formula syntax | Method |
|----|----|----|
| Continuous measurement error | `me(x, sd)` | SIMEX (Cook & Stefanski, 1994) |
| Berkson measurement error | `me(x, sd, type = "berkson")` | Berkson SIMEX |
| Systematic (non-zero mean) ME | `me(x, sd, mean = mu)` | Non-zero mean SIMEX |
| Discrete misclassification | `mc(z, Pi)` | MC-SIMEX (Kuchenhoff et al., 2006) |
| Improved misclassification | `mc(z, Pi)` with `method = "improved"` | Sevilimedu & Yu (2026) |
| Response misclassification | `mc(y, Pi_y) ~ ...` | MC-SIMEX on response |

### Supported families

Gaussian ([`gaussian()`](https://rdrr.io/r/stats/family.html)), Poisson
([`poisson()`](https://rdrr.io/r/stats/family.html)), and Binomial
([`binomial()`](https://rdrr.io/r/stats/family.html)).

## Background: the SIMEX algorithm

When a covariate $`X`$ is measured with additive error $`W = X + U`$,
where $`U \sim N(0, \sigma_u^2)`$, regression coefficients are biased
toward zero (attenuation bias). SIMEX corrects this in two steps:

1.  **Simulation**: For $`\lambda \in \{0.5, 1, 1.5, 2\}`$ and
    $`b = 1, \ldots, B`$, generate
    $`W^*_{b,\lambda} = W + \sqrt{\lambda} \cdot U_b`$, where
    $`U_b \sim N(0, \sigma_u^2)`$. Refit the model on each
    pseudo-dataset.

2.  **Extrapolation**: Fit a polynomial (quadratic, linear, or
    loglinear) to the average estimates $`\bar{\hat\beta}(\lambda)`$ as
    a function of $`\lambda`$, and extrapolate to $`\lambda = -1`$ (no
    measurement error).

For misclassification, the MC-SIMEX algorithm replaces additive noise
with resampling from $`\Pi^\lambda`$, the misclassification matrix
raised to fractional powers.

## Example 1: Continuous measurement error

``` r

n <- 2000
x_true <- rnorm(n)
y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
x_obs <- x_true + rnorm(n, sd = 0.5)  # measurement error with sd = 0.5
df <- data.frame(y = y, x = x_obs)

# Naive OLS (biased)
naive <- lm(y ~ x, data = df)

# SIMEX corrected
fit <- simex(y ~ me(x, 0.5), data = df, B = 200, seed = 42)

cat("True slope:  2.000\n")
#> True slope:  2.000
cat("Naive slope:", round(coef(naive)["x"], 3), "\n")
#> Naive slope: 1.578
cat("SIMEX slope:", round(coef(fit)["x"], 3), "\n")
#> SIMEX slope: 1.915
```

The naive estimate is attenuated toward zero. SIMEX recovers a value
closer to the true slope of 2.

### Summary and inference

``` r

summary(fit)
#> 
#> Call:
#> simex(formula = y ~ me(x, 0.5), data = df, B = 200, seed = 42)
#> 
#> Family: gaussian 
#> SIMEX variable(s): x 
#> Extrapolation: quadratic 
#> Lambda grid: 0, 0.5, 1, 1.5, 2 
#> B = 200 , n = 2000 
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -3.38970 -0.78080 -0.02356  0.79428  3.25930 
#> 
#> Naive coefficients:
#> (Intercept)           x 
#>      0.9843      1.5779 
#> 
#> SIMEX corrected coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.98948    0.02418   40.92   <2e-16 ***
#> x            1.91518    0.02453   78.08   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
confint(fit)
#>                 2.5 %   97.5 %
#> (Intercept) 0.9420865 1.036866
#> x           1.8671012 1.963252
```

### SIMEX trajectory plot

``` r

plot(fit)
```

![](mismeasured_files/figure-html/me-plot-1.png)![](mismeasured_files/figure-html/me-plot-2.png)

The plot shows coefficient estimates at each $`\lambda`$ level (solid
points), the naive estimate at $`\lambda = 0`$ (open circle), the
extrapolation curve, and the corrected estimate at $`\lambda = -1`$ (red
cross).

### Changing the extrapolation method

``` r

fit_lin <- refit(fit, extrapolation = "linear")
cat("Quadratic:", round(coef(fit)["x"], 3), "\n")
#> Quadratic: 1.915
cat("Linear:   ", round(coef(fit_lin)["x"], 3), "\n")
#> Linear:    1.78
```

## Example 2: Multiple error-prone covariates

``` r

x1_true <- rnorm(n)
x2_true <- rnorm(n)
y <- 1 + 1.5 * x1_true - 0.8 * x2_true + rnorm(n, sd = 0.5)
df2 <- data.frame(
  y  = y,
  x1 = x1_true + rnorm(n, sd = 0.3),
  x2 = x2_true + rnorm(n, sd = 0.5)
)

fit2 <- simex(y ~ me(x1, 0.3) + me(x2, 0.5), data = df2, B = 200, seed = 42)
summary(fit2)
#> 
#> Call:
#> simex(formula = y ~ me(x1, 0.3) + me(x2, 0.5), data = df2, B = 200, 
#>     seed = 42)
#> 
#> Family: gaussian 
#> SIMEX variable(s): x1, x2 
#> Extrapolation: quadratic 
#> Lambda grid: 0, 0.5, 1, 1.5, 2 
#> B = 200 , n = 2000 
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -2.431169 -0.525458  0.008814  0.508883  2.844247 
#> 
#> Naive coefficients:
#> (Intercept)          x1          x2 
#>      0.9901      1.3756     -0.6430 
#> 
#> SIMEX corrected coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.99644    0.01751   56.92   <2e-16 ***
#> x1           1.48999    0.01691   88.13   <2e-16 ***
#> x2          -0.77492    0.01739  -44.55   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Example 3: Misclassification (standard MC-SIMEX)

A binary covariate $`Z`$ is misclassified with known misclassification
matrix $`\Pi`$, where $`\Pi_{jl} = P(Z^* = j \mid Z = l)`$.

``` r

z <- rbinom(n, 1, 0.4)
y3 <- 1 + 2 * z + 0.5 * rnorm(n) + rnorm(n)

# Misclassify: sensitivity = 0.85, specificity = 0.9
z_star <- z
z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
df3 <- data.frame(y = y3, z = factor(z_star))

fit3 <- simex(y ~ mc(z, Pi), data = df3, method = "standard",
              B = 200, seed = 42)

cat("True z effect:     2.000\n")
#> True z effect:     2.000
cat("Naive z effect:   ", round(fit3$naive.coefficients[1], 3), "\n")
#> Naive z effect:    1.485
cat("MC-SIMEX z effect:", round(coef(fit3)[1], 3), "\n")
#> MC-SIMEX z effect: 1.962
```

## Example 4: Improved MC-SIMEX (exact correction)

The improved method (Sevilimedu & Yu, 2026) computes an exact
fixed-matrix correction instead of extrapolation. It requires only B = 1
replicate and is the default for a single misclassified covariate. For
binary covariates, the scalar correction factor is reported as
`c.lambda`; for K-level factors, the dummy-coefficient correction matrix
is stored in `correction.matrix`.

``` r

fit4 <- simex(y ~ mc(z, Pi), data = df3, seed = 42)  # default: method = "improved"

cat("Method:", fit4$method, "\n")
#> Method: improved
cat("Correction factor c(lambda):", round(fit4$c.lambda, 4), "\n")
#> Correction factor c(lambda): 1.7954
cat("Improved z effect:", round(coef(fit4)[1], 3), "\n")
#> Improved z effect: 2.001
```

### Optimal lambda

The improved method can search for the lambda that minimizes the
absolute correction factor:

``` r

fit4_opt <- simex(y ~ mc(z, Pi), data = df3, lambda = "optimal", seed = 42)
cat("Optimal lambda:", fit4_opt$lambda[-1], "\n")
#> Optimal lambda: 0.1
cat("Optimal c(lambda):", round(fit4_opt$c.lambda, 4), "\n")
#> Optimal c(lambda): 1.3808
```

## Example 5: Multiple misclassified covariates

``` r

z1 <- rbinom(n, 1, 0.4)
z2 <- rbinom(n, 1, 0.6)
y5 <- 1 + 2 * z1 + 1.5 * z2 + rnorm(n)

# Misclassify both
z1_star <- z1
z1_star[z1 == 0] <- rbinom(sum(z1 == 0), 1, 0.10)
z1_star[z1 == 1] <- 1 - rbinom(sum(z1 == 1), 1, 0.15)

z2_star <- z2
z2_star[z2 == 0] <- rbinom(sum(z2 == 0), 1, 0.08)
z2_star[z2 == 1] <- 1 - rbinom(sum(z2 == 1), 1, 0.12)

Pi1 <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
Pi2 <- matrix(c(0.92, 0.08, 0.12, 0.88), 2, 2)
df5 <- data.frame(y = y5, z1 = factor(z1_star), z2 = factor(z2_star))

fit5 <- simex(y ~ mc(z1, Pi1) + mc(z2, Pi2), data = df5,
              method = "standard", B = 200, seed = 42)
summary(fit5)
#> 
#> Call:
#> simex(formula = y ~ mc(z1, Pi1) + mc(z2, Pi2), data = df5, method = "standard", 
#>     B = 200, seed = 42)
#> 
#> Family: gaussian 
#> MC-SIMEX variables: z1, z2 
#> Method: standard 
#> Extrapolation: quadratic 
#> Lambda grid: 0, 0.5, 1, 1.5, 2 
#> B = 200 , n = 2000 
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -4.421865 -0.789852  0.005384  0.824377  4.935690 
#> 
#> Naive coefficients:
#>         z11         z21 (Intercept) 
#>      1.4812      1.2013      1.4261 
#> 
#> MC-SIMEX corrected coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> z11          1.91470    0.08737   21.91   <2e-16 ***
#> z21          1.47276    0.08430   17.47   <2e-16 ***
#> (Intercept)  1.08076    0.06494   16.64   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

With multiple
[`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) terms,
the standard method is used automatically. Each variable is resampled
independently via its own misclassification matrix.

## Example 6: Response misclassification

When the outcome itself is misclassified (e.g., disease status from an
imperfect diagnostic test), use
[`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) on the
left-hand side:

``` r

x6 <- rnorm(n)
y6_true <- rbinom(n, 1, plogis(-0.5 + 1.2 * x6))

# Misclassify the response
y6_obs <- y6_true
y6_obs[y6_true == 0] <- rbinom(sum(y6_true == 0), 1, 0.10)
y6_obs[y6_true == 1] <- 1 - rbinom(sum(y6_true == 1), 1, 0.08)

Pi_y <- matrix(c(0.90, 0.10, 0.08, 0.92), 2, 2)
df6 <- data.frame(y = factor(y6_obs), x = x6)

fit6 <- simex(mc(y, Pi_y) ~ x, family = binomial(), data = df6,
              B = 200, seed = 42)

cat("True x effect:     1.200\n")
#> True x effect:     1.200
cat("Naive x effect:   ", round(fit6$naive.coefficients["x"], 3), "\n")
#> Naive x effect:    0.9
cat("MC-SIMEX x effect:", round(coef(fit6)["x"], 3), "\n")
#> MC-SIMEX x effect: 1.175
```

Response misclassification can be combined with covariate
misclassification:

``` r

z6 <- rbinom(n, 1, 0.4)
y6b_true <- rbinom(n, 1, plogis(-0.5 + 1.0 * z6 + 0.5 * x6))
y6b_obs <- y6b_true
y6b_obs[y6b_true == 0] <- rbinom(sum(y6b_true == 0), 1, 0.10)
y6b_obs[y6b_true == 1] <- 1 - rbinom(sum(y6b_true == 1), 1, 0.08)

z6_star <- z6
z6_star[z6 == 0] <- rbinom(sum(z6 == 0), 1, 0.10)
z6_star[z6 == 1] <- 1 - rbinom(sum(z6 == 1), 1, 0.15)

df6b <- data.frame(y = factor(y6b_obs), z = factor(z6_star), x = x6)

fit6b <- simex(mc(y, Pi_y) ~ mc(z, Pi) + x, family = binomial(),
               data = df6b, B = 200, seed = 42)
summary(fit6b)
#> 
#> Call:
#> simex(formula = mc(y, Pi_y) ~ mc(z, Pi) + x, family = binomial(), 
#>     data = df6b, B = 200, seed = 42)
#> 
#> Family: binomial 
#> MC-SIMEX variable: z 
#> Method: standard 
#> Extrapolation: quadratic 
#> Lambda grid: 0, 0.5, 1, 1.5, 2 
#> B = 200 , n = 2000 
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -0.8432 -0.4503 -0.2201  0.4895  0.8450 
#> 
#> Naive coefficients:
#>          z1 (Intercept)           x 
#>      0.4873     -0.2085      0.3293 
#> 
#> MC-SIMEX corrected coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> z1           0.75905    0.13445   5.646 1.88e-08 ***
#> (Intercept) -0.35817    0.07936  -4.513 6.76e-06 ***
#> x            0.39681    0.05843   6.791 1.46e-11 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Example 7: Berkson measurement error

In the Berkson error model, the true value is a noisy version of the
observed (assigned) value: $`X = W + U`$. This arises in occupational
exposure assessment, where $`W`$ is the assigned group mean and $`X`$ is
the individual’s true exposure.

``` r

w <- rnorm(n)
u <- rnorm(n, sd = 0.5)
x_true_b <- w + u
y7 <- 1 + 2 * x_true_b + rnorm(n, sd = 0.5)
df7 <- data.frame(y = y7, x = w)

fit7_class <- simex(y ~ me(x, 0.5), data = df7, B = 200, seed = 42)
fit7_berk  <- simex(y ~ me(x, 0.5, type = "berkson"), data = df7,
                    B = 200, seed = 42)

cat("Classical SIMEX slope:", round(coef(fit7_class)["x"], 3), "\n")
#> Classical SIMEX slope: 2.567
cat("Berkson SIMEX slope:  ", round(coef(fit7_berk)["x"], 3), "\n")
#> Berkson SIMEX slope:   2.562
```

## Example 8: Non-zero mean measurement error

When measurement error has a systematic component
$`U \sim N(\mu, \sigma^2)`$ with $`\mu \neq 0`$:

``` r

x_true_m <- rnorm(n)
x_obs_m <- x_true_m + rnorm(n, mean = 0.3, sd = 0.5)
y8 <- 1 + 2 * x_true_m + rnorm(n, sd = 0.5)
df8 <- data.frame(y = y8, x = x_obs_m)

fit8_std <- simex(y ~ me(x, 0.5), data = df8, B = 200, seed = 42)
fit8_nzm <- simex(y ~ me(x, 0.5, mean = 0.3), data = df8, B = 200, seed = 42)

cat("Standard SIMEX (ignoring mean):", round(coef(fit8_std)["x"], 3), "\n")
#> Standard SIMEX (ignoring mean): 1.938
cat("Non-zero mean SIMEX:           ", round(coef(fit8_nzm)["x"], 3), "\n")
#> Non-zero mean SIMEX:            1.938
```

## Example 9: Poisson regression with misclassification

``` r

z9 <- rbinom(n, 1, 0.4)
y9 <- rpois(n, exp(0.5 + 0.8 * z9))

z9_star <- z9
z9_star[z9 == 0] <- rbinom(sum(z9 == 0), 1, 0.10)
z9_star[z9 == 1] <- 1 - rbinom(sum(z9 == 1), 1, 0.15)

df9 <- data.frame(y = y9, z = factor(z9_star))

fit9 <- simex(y ~ mc(z, Pi), family = poisson(), data = df9, seed = 42)
summary(fit9)
#> 
#> Call:
#> simex(formula = y ~ mc(z, Pi), family = poisson(), data = df9, 
#>     seed = 42)
#> 
#> Family: poisson 
#> MC-SIMEX variable: z 
#> Method: improved 
#> Extrapolation: exact (improved) 
#> Lambda grid: 0, 1 
#> B = 1 , n = 2000 
#> Estimated P(X=1): 0.3807 
#> Correction factor(s): 1.7925 
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -3.8076 -1.7772  0.1924  1.1924  7.2228 
#> 
#> Naive coefficients:
#>           1 (Intercept) 
#>      0.5784      0.6237 
#> 
#> MC-SIMEX corrected coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> 1            0.76194    0.05151   14.79   <2e-16 ***
#> (Intercept)  0.57505    0.02762   20.82   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Example 10: Frequency weights

The `weights` argument supports frequency (prior) weights:

``` r

# Aggregated data with weights
df10_agg <- data.frame(
  y = c(2.5, 3.1, 4.2, 5.0),
  x = c(1.0, 2.0, 3.0, 4.0),
  freq = c(100, 150, 120, 130)
)

fit10 <- simex(y ~ me(x, 0.3), data = df10_agg, weights = df10_agg$freq,
               B = 100, seed = 42)
coef(fit10)
#> (Intercept)           x 
#>   1.4042179   0.9413936
```

## Utility functions

### Building misclassification matrices

[`build.mc.matrix()`](https://ojalab.github.io/mismeasured/reference/build.mc.matrix.md)
finds the nearest valid misclassification matrix (one that can be raised
to fractional powers without negative probabilities):

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

### Checking matrix validity

``` r

P1 <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2)
P2 <- matrix(c(0.4, 0.6, 0.6, 0.4), 2, 2)  # too much error
check.mc.matrix(list(P1, P2))
#> [1]  TRUE FALSE
```

## Implementation details

### C++ backend

The simulation step runs entirely in C++ for all error types:

- `simex_sim_cpp`: Continuous ME with IRLS fitting + model vcov
  accumulation
- `mcsimex_sim_cpp`: Single-variable MC-SIMEX
- `mcsimex_multi_sim_cpp`: Multi-variable MC-SIMEX with optional
  response mc

The IRLS solver supports Gaussian (identity link), Poisson (log link),
and Binomial (logit link) families. Matrix exponentiation for
$`\Pi^\lambda`$ uses eigendecomposition.

### Variance estimation

- **ME-SIMEX**: Jackknife variance following Stefanski & Cook (1994).
  The model-based variance $`\bar{V}_{\text{model}}(\lambda)`$ is the
  average of $`(X_{\text{sim}}'WX_{\text{sim}})^{-1}\phi`$ across B
  refits, computed in C++.
  $`V_{\text{jk}}(\lambda) = \bar{V}_{\text{model}} - \text{Cov}(\hat\theta_b)`$
  is then extrapolated to $`\lambda = -1`$.

- **MC-SIMEX standard**: Same jackknife approach with sandwich variance
  at the mean estimate per lambda.

- **MC-SIMEX improved**: Scaled model-based variance:
  $`V_{\text{corrected}} = c(\lambda)^2 \cdot V_{\text{model}}`$.

## References

Cook, J.R. and Stefanski, L.A. (1994). Simulation-extrapolation
estimation in parametric measurement error models. *Journal of the
American Statistical Association*, 89(428), 1314-1328.

Kuchenhoff, H., Mwalili, S.M. and Lesaffre, E. (2006). A general method
for dealing with misclassification in regression: The misclassification
SIMEX. *Biometrics*, 62(1), 85-96.

Sevilimedu, V. and Yu, L. (2026). An improved misclassification
simulation extrapolation (MC-SIMEX) algorithm. *Statistics in Medicine*,
45, e70418.
