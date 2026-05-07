# SIMEX estimator for GLMs with measurement error or misclassification

Fits a generalized linear model corrected for measurement error or
misclassification using the Simulation-Extrapolation (SIMEX) method.
Error-prone variables are specified directly in the formula using
[`me`](https://ojalab.github.io/mismeasured/reference/me.md) (continuous
measurement error) or
[`mc`](https://ojalab.github.io/mismeasured/reference/mc.md) (discrete
misclassification).

## Usage

``` r
simex(
  formula,
  family = gaussian(),
  data,
  method = NULL,
  lambda = NULL,
  B = NULL,
  extrapolation = c("quadratic", "linear", "loglinear"),
  jackknife = TRUE,
  weights = NULL,
  seed = 42L
)
```

## Arguments

- formula:

  a formula with
  [`me`](https://ojalab.github.io/mismeasured/reference/me.md) and/or
  [`mc`](https://ojalab.github.io/mismeasured/reference/mc.md) terms.
  See Examples.

- family:

  a GLM family (default:
  [`gaussian()`](https://rdrr.io/r/stats/family.html)).

- data:

  a data frame.

- method:

  character: correction method for
  [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) terms.
  `NULL` (default) auto-selects `"improved"`. Ignored for
  [`me()`](https://ojalab.github.io/mismeasured/reference/me.md) terms.
  See Details.

- lambda:

  numeric vector of SIMEX exponents, or `"optimal"` (improved MC-SIMEX
  only). Default: auto-set based on error type.

- B:

  integer number of simulation replicates. Default: auto-set.

- extrapolation:

  extrapolation method: `"quadratic"` (default), `"linear"`, or
  `"loglinear"`. Ignored for improved MC-SIMEX.

- jackknife:

  logical: compute variance estimates? Default `TRUE`.

- weights:

  optional prior weights.

- seed:

  random seed (default: 42).

## Value

An object of class `"simex"`.

## Details

The function auto-detects the error type from the formula:

- **Continuous measurement error**
  ([`me()`](https://ojalab.github.io/mismeasured/reference/me.md)
  terms): Uses the SIMEX algorithm of Cook and Stefanski (1994). Extra
  noise is added to the error-prone variable, the model is refitted B
  times per lambda level, and coefficients are extrapolated to lambda =
  -1.

- **Discrete misclassification**
  ([`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md)
  terms): Uses the MC-SIMEX algorithm. With `method = "improved"`
  (default), the exact fixed-matrix correction of Sevilimedu and
  Yu (2026) and its K-level dummy-vector extension are applied for one
  misclassified covariate, requiring only B = 1 replicate. With
  `method = "standard"`, the original Kuchenhoff et al. (2006)
  extrapolation-based approach is used.

## References

Cook, J.R. and Stefanski, L.A. (1994). Simulation-extrapolation
estimation in parametric measurement error models. *JASA*, 89,
1314–1328.

Kuchenhoff, H., Mwalili, S.M. and Lesaffre, E. (2006). A general method
for dealing with misclassification in regression: The misclassification
SIMEX. *Biometrics*, 62(1), 85–96.

Sevilimedu, V. and Yu, L. (2026). An improved misclassification
simulation extrapolation (MC-SIMEX) algorithm. *Statistics in Medicine*,
45, e70418.

## Examples

``` r
# --- Continuous measurement error ---
set.seed(42)
n <- 500
x_true <- rnorm(n)
y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
x_obs <- x_true + rnorm(n, sd = 0.5)
df <- data.frame(y = y, x = x_obs)

fit_me <- simex(y ~ me(x, 0.5), data = df, B = 50)
summary(fit_me)
#> 
#> Call:
#> simex(formula = y ~ me(x, 0.5), data = df, B = 50)
#> 
#> Family: gaussian 
#> SIMEX variable(s): x 
#> Extrapolation: quadratic 
#> Lambda grid: 0, 0.5, 1, 1.5, 2 
#> B = 50 , n = 500 
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -3.481486 -0.731792 -0.005704  0.676633  3.281053 
#> 
#> Naive coefficients:
#> (Intercept)           x 
#>      1.0094      1.5814 
#> 
#> SIMEX corrected coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  1.02100    0.03932   25.96   <2e-16 ***
#> x            1.93753    0.04507   42.99   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# --- Misclassification (improved, default) ---
Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
z <- rbinom(n, 1, 0.4)
y2 <- rpois(n, exp(0.5 + 0.8 * z + 0.3 * x_true))
z_star <- z
z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
df2 <- data.frame(y = y2, z = factor(z_star), x = x_true)

fit_mc <- simex(y ~ mc(z, Pi) + x, family = poisson(), data = df2)
summary(fit_mc)
#> 
#> Call:
#> simex(formula = y ~ mc(z, Pi) + x, family = poisson(), data = df2)
#> 
#> Family: poisson 
#> MC-SIMEX variable: z 
#> Method: improved 
#> Extrapolation: exact (improved) 
#> Lambda grid: 0, 1 
#> B = 1 , n = 500 
#> Estimated P(X=1): 0.416 
#> Correction factor(s): 1.7688 
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -4.7607 -1.3396 -0.2456  0.9407  6.3076 
#> 
#> Naive coefficients:
#>           1 (Intercept)           x 
#>      0.5258      0.6865      0.2914 
#> 
#> MC-SIMEX corrected coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> 1            0.74272    0.09785   7.591 1.59e-13 ***
#> (Intercept)  0.59943    0.05581  10.740  < 2e-16 ***
#> x            0.29214    0.02804  10.418  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
```
