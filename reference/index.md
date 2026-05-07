# Package index

## Bias-corrected GLM / LM

Analytical bias-correction estimators (BCA, BCM, CS) and the joint
mixture-likelihood one-step estimator for GLMs with a misclassified
covariate.

- [`mcglm()`](https://ojalab.github.io/mismeasured/reference/mcglm.md) :
  Bias-Corrected GLM Estimation with Misclassified Covariates
- [`mclm()`](https://ojalab.github.io/mismeasured/reference/mclm.md) :
  Bias-Corrected Linear Model with Misclassified Covariates

## SIMEX and MC-SIMEX

Simulation-extrapolation for measurement error and misclassification.

- [`simex()`](https://ojalab.github.io/mismeasured/reference/simex.md) :
  SIMEX estimator for GLMs with measurement error or misclassification
- [`misclass()`](https://ojalab.github.io/mismeasured/reference/misclass.md)
  : Generate misclassified data
- [`me()`](https://ojalab.github.io/mismeasured/reference/me.md) :
  Specify a variable measured with error
- [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) :
  Specify a misclassified variable

## Misclassification-matrix utilities

- [`build.mc.matrix()`](https://ojalab.github.io/mismeasured/reference/build.mc.matrix.md)
  : Build a valid misclassification matrix
- [`check.mc.matrix()`](https://ojalab.github.io/mismeasured/reference/check.mc.matrix.md)
  : Check if a misclassification matrix is valid for fractional powers
- [`diag.block()`](https://ojalab.github.io/mismeasured/reference/diag.block.md)
  : Construct a block diagonal matrix

## Methods for mcglm fits

S3 methods for objects returned by mcglm() and mclm().

- [`summary(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/summary.mcglm.md)
  :

  Summary of an `mcglm` fit

- [`print(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/print.mcglm.md)
  :

  Print method for `mcglm` fits

- [`coef(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/coef.mcglm.md)
  :

  Extract coefficients from an `mcglm` fit

- [`vcov(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/vcov.mcglm.md)
  :

  Variance-covariance matrix of an `mcglm` fit

- [`confint(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/confint.mcglm.md)
  :

  Wald confidence intervals for an `mcglm` fit

- [`se.mcglm()`](https://ojalab.github.io/mismeasured/reference/se.mcglm.md)
  : Standard errors of coefficient estimates

- [`predict(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/predict.mcglm.md)
  : Linear-predictor or response predictions

- [`fitted(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/fitted.mcglm.md)
  :

  Fitted values for an `mcglm` fit

- [`residuals(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/residuals.mcglm.md)
  :

  Residuals for an `mcglm` fit

- [`logLik(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/logLik.mcglm.md)
  :

  Log-likelihood of an `mcglm` fit

- [`AIC(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/AIC.mcglm.md)
  :

  AIC of an `mcglm` fit

- [`BIC(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/BIC.mcglm.md)
  :

  BIC of an `mcglm` fit

- [`nobs(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/nobs.mcglm.md)
  : Number of observations

- [`family(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/family.mcglm.md)
  :

  Family of an `mcglm` fit

- [`formula(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/formula.mcglm.md)
  :

  Formula of an `mcglm` fit (NULL for matrix interface)

- [`model.matrix(`*`<mcglm>`*`)`](https://ojalab.github.io/mismeasured/reference/model.matrix.mcglm.md)
  :

  Proxy design matrix used by `mcglm`

## Other methods

- [`refit()`](https://ojalab.github.io/mismeasured/reference/refit.md) :
  Refit a SIMEX model with a different extrapolation function
