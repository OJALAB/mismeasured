# nonprobsvyMC 0.1.0

* Initial version: `nonprob_mc()`, a mass-imputation estimator of the
  population mean with a misclassification-corrected outcome model
  (`mismeasured::mcglm()`; naive/BCA/BCM/CS via `mc_method`), returning
  `nonprob`-class objects compatible with the `nonprobsvy` methods
  (`print`, `summary`, `confint`, `nobs`, `weights`, `pop_size`).
* Analytic variance: design-based `var_prob` from **survey** plus the
  corrected-score linearization `var_nonprob` assembled from
  `sandwich::bread()`/`meat()` of the outcome fit.
* Supports binary and multicategory (K >= 3) misclassified covariates
  with known misclassification parameters.
