# Variance estimation for fixed-Pi K-level improved MC-SIMEX

Sampling variance: average over lambda of T_l V(theta(lambda_l)) T_l',
where V(theta(lambda_l)) is the vcov of a refit on a Pi^lambda_l
resampled design (Sevilimedu & Yu 2026, eq. 4, generalised to the
K-level transform). For n_lambda \* B \> 1 a between-replicate
Monte-Carlo term is added; cov(all_corrected) alone is only the
resampling variability of the mean and would shrink to zero as B grows.

## Usage

``` r
.variance_k_improved(theta_list, vcov_list, transform_list, lambda, B, p)
```
