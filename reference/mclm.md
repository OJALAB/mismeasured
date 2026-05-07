# Bias-Corrected Linear Model with Misclassified Covariates

Convenience wrapper around
[`mcglm`](https://ojalab.github.io/mismeasured/reference/mcglm.md) that
fixes the family to `"gaussian"` (identity link), giving a
misclassification-corrected linear-regression interface analogous to
[`lm`](https://rdrr.io/r/stats/lm.html). All arguments other than
`family` are forwarded to `mcglm` and carry the same meaning.

## Usage

``` r
mclm(
  formula,
  data = NULL,
  method = c("naive", "bca", "bcm", "cs"),
  p01 = NULL,
  p10 = NULL,
  pi_z = NULL,
  Pi = NULL,
  K = NULL,
  c1 = NULL,
  c2 = NULL,
  iterate = FALSE,
  jacobian = c("analytical", "numerical"),
  fix_omega = FALSE,
  vcov_corrected = FALSE,
  weights = NULL,
  J = NULL,
  homoskedastic = TRUE,
  optim_control = list(),
  z_hat = NULL,
  x = NULL
)
```

## Arguments

- formula:

  Either a formula of the form `y ~ mc(z, Pi) + x1 + x2`, where
  [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) marks
  the misclassified covariate \\\hat Z_i\\ together with its
  misclassification matrix \\\Pi\\; or, when using the matrix interface,
  a numeric response vector of length \\n\\.

- data:

  A data frame (required for the formula interface). The variable named
  in [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md)
  should be a [`factor`](https://rdrr.io/r/base/factor.html) or
  integer-valued vector with levels coded as \\0, 1, \dots, K-1\\.

- method:

  Character vector of estimators to fit. Any subset of
  `c("naive", "bca", "bcm", "cs", "cs_akn", "onestep")`; the default is
  the four analytical estimators `c("naive", "bca", "bcm", "cs")`. The
  `"cs_akn"` entry selects the Akazawa–Kinukawa–Nakamura (1998)
  corrected-score construction (an alternative formulation of `"cs"`
  based on the unbiased-surrogate transform \\x = Q^{-1}(u - p_0)\\); it
  needs only \\\Pi\\ (not \\\pi_z\\) and is unsupported for
  `family = "multinomial"`.

- p01:

  False-positive rate \\p\_{01} = \Pr(\hat Z = 1 \mid Z = 0)\\ (`K = 2`
  only). Auto-extracted from `Pi` when supplied.

- p10:

  False-negative rate \\p\_{10} = \Pr(\hat Z = 0 \mid Z = 1)\\ (`K = 2`
  only). Auto-extracted from `Pi` when supplied.

- pi_z:

  Latent prevalence: a scalar \\\pi = \Pr(Z = 1)\\ (`K = 2`) or a
  \\K\\-vector \\(\pi_0, \dots, \pi\_{K-1})\\ of class probabilities. If
  `NULL`, `pi_z` is estimated from `z_hat` and `Pi` via \\\hat\pi =
  \Pi^{-1} \hat\pi\_{\text{obs}}\\ (clamped to \\\[0.01, 0.99\]\\).

- Pi:

  The \\K \times K\\ misclassification matrix \\\Pi\_{j\ell} = \Pr(\hat
  Z = j - 1 \mid Z = \ell - 1)\\; columns must sum to 1. Extracted
  automatically from the
  [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) term
  when using the formula interface.

- K:

  Number of latent categories. Auto-detected from `Pi` or, failing that,
  from the levels of `z_hat`.

- c1:

  Pre-computed misclassification constant \\c_1 = p\_{01}(1 - \pi)\\
  (binary case, advanced use). Computed from `p01` and `pi_z` when
  `NULL`.

- c2:

  Pre-computed misclassification constant \\c_2 = p\_{01}(1 - \pi) -
  p\_{10} \pi\\ (binary case, advanced use). Computed from `p01`, `p10`
  and `pi_z` when `NULL`.

- iterate:

  Logical. If `TRUE`, iterate the BCA/BCM updates to convergence (the
  iterated BCM solves the corrected score equation; iterated BCA
  performs Newton steps on the additive correction).

- jacobian:

  Character. Either `"analytical"` (default) or `"numerical"`. Controls
  how the drift Jacobian \\M = \partial m / \partial \psi\\ is evaluated
  for the multicategory (\\K \ge 3\\) BCM and CS estimators. The
  analytical form is exact and faster; the numerical form (forward
  differences) is provided as a fallback for cross-checking. The binary
  case (\\K = 2\\) always uses the analytical Jacobian.

- fix_omega:

  Logical. If `TRUE`, fix the mixture weights in the one-step estimator
  at the values implied by the supplied misclassification parameters; if
  `FALSE` (default), they are estimated jointly with \\\psi\\ via a
  softmax parameterization.

- vcov_corrected:

  Logical. If `TRUE`, the BCA/BCM variance accounts for the additional
  uncertainty due to plug-in estimation of \\\hat m(\hat\psi)\\ (joint
  score-and-drift sandwich); if `FALSE` (default), the drifting-regime
  asymptotic variance \\A^{-1} C A^{-1}\\ is used. Both share the same
  point estimate.

- weights:

  Optional positive frequency weights of length \\n\\.

- J:

  Number of response categories (multinomial family only;
  auto-detected).

- homoskedastic:

  Logical. For one-step Gaussian fits, assume a single residual variance
  across mixture components; if `FALSE`, separate variances are
  estimated.

- optim_control:

  List of control parameters passed to
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html) for the one-step
  estimator.

- z_hat:

  Integer vector of observed proxy values \\\hat Z_i \in \\0, \dots,
  K-1\\\\ (matrix interface only).

- x:

  Covariate matrix \\\[x_1, \dots, x_n\]^\top\\ of dimension \\n \times
  r\\; should include an intercept column when one is desired (matrix
  interface only).

## Value

An object of class `"mcglm"` with components:

- coefficients:

  Named list of coefficient vectors \\\hat\psi\\, one per method.

- vcov:

  Named list of asymptotic variance-covariance matrices, one per method
  (sandwich estimators from the paper's theorems).

- se:

  Named list of standard errors per method.

- naive_fit:

  The [`glm`](https://rdrr.io/r/stats/glm.html) object from the naive
  fit (`NULL` for multinomial).

- method, family, K, n, p:

  Metadata.

- weights:

  Frequency weights used (`NULL` if unweighted).

- loglik_onestep, vcov_onestep:

  One-step log-likelihood and variance (when `onestep` is fit).

- call, formula:

  The matched call and the formula (`NULL` if the matrix interface was
  used).

See
[`summary.mcglm`](https://ojalab.github.io/mismeasured/reference/summary.mcglm.md),
[`vcov.mcglm`](https://ojalab.github.io/mismeasured/reference/vcov.mcglm.md),
[`confint.mcglm`](https://ojalab.github.io/mismeasured/reference/confint.mcglm.md),
[`predict.mcglm`](https://ojalab.github.io/mismeasured/reference/predict.mcglm.md),
and
[`logLik.mcglm`](https://ojalab.github.io/mismeasured/reference/logLik.mcglm.md)
for downstream methods.

## See also

[`mcglm`](https://ojalab.github.io/mismeasured/reference/mcglm.md) for
the general GLM interface and the underlying methodology;
[`simex`](https://ojalab.github.io/mismeasured/reference/simex.md) for
SIMEX/MC-SIMEX alternatives.

## Examples

``` r
set.seed(33)
n  <- 1000
Pi <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
z  <- rbinom(n, 1, 0.4)
z_hat <- z
z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.1)
x1 <- rnorm(n)
y  <- -0.5 + 0.8 * z + 0.7 * x1 + rnorm(n, 0, 0.5)
df <- data.frame(y = y, z = factor(z_hat), x1 = x1)

fit <- mclm(y ~ mc(z, Pi) + x1, data = df)
summary(fit)
#> 
#> Call:
#> mclm(formula = y ~ mc(z, Pi) + x1, data = df)
#> 
#> Family: gaussian  |  n = 1000, K = 2, p = 3
#> Methods: naive, bca, bcm, cs
#> 
#> --- NAIVE ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.61931    0.03542   17.48   <2e-16 ***
#> alpha0 -0.44109    0.02203  -20.03   <2e-16 ***
#> alpha1  0.69747    0.01758   39.66   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> --- BCA ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.75130    0.03574   21.02   <2e-16 ***
#> alpha0 -0.48301    0.02210  -21.86   <2e-16 ***
#> alpha1  0.69865    0.01780   39.25   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> --- BCM ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.78706    0.03593   21.91   <2e-16 ***
#> alpha0 -0.49436    0.02214  -22.33   <2e-16 ***
#> alpha1  0.69897    0.01790   39.05   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> --- CS ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.78706    0.04516   17.43   <2e-16 ***
#> alpha0 -0.49436    0.02407  -20.54   <2e-16 ***
#> alpha1  0.69897    0.01786   39.13   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual deviance (naive): 299.1  on 997 degrees of freedom
#> AIC (naive): 1639
#> 
#> Bias correction (difference from naive):
#>         bca        bcm        cs       
#> gamma    0.131998   0.167752   0.167752
#> alpha0  -0.041921  -0.053276  -0.053276
#> alpha1   0.001178   0.001498   0.001498
```
