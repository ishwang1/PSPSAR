# PSPSAR

`PSPSAR` estimates spatial autoregressive models in which the spatial
spillover coefficient can vary with a common observed state.  It also provides
a constant coefficient SAR benchmark.

## Installation

Install the public GitHub version with `remotes`:

```r
install.packages("remotes")
remotes::install_github("ishwang1/PSPSAR")
```

The package requires R (>= 4.1.0), `MASS`, and `sandwich`.

## Models

- `PSPSAR()` fits a profiled state-varying SAR model. It supports 2SLS or
  feasible GMM estimation, debiased or undersmoothed inference, several kernel
  choices, and Newey-West or Andrews HAC covariance estimation.
- `SAR()` fits the constant coefficient SAR benchmark using 2SLS or feasible
  GMM with HAC inference.

## Required inputs

Users supply their own data. `Y` must be a numeric `T` by `N` outcome matrix.
`X` may be an `N` by `p` time-invariant covariate matrix or a `T` by `N` by `p`
array of time-varying covariates. `W` must be an `N` by `N` spatial weight
matrix. `PSPSAR()` additionally requires a numeric state vector `z` of length
`T`.

Both functions add an intercept internally. Set `preprocess_W = TRUE` only
when `W` has already been prepared for estimation; otherwise the functions
apply their built-in spatial weight preprocessing.

Spatial-lag instruments are generated only for covariates that exhibit
cross-sectional variation in at least one period. A covariate that is common
to all units within every period, such as a time fixed effect, remains in `X` but does not receive a redundant `W X`
instrument. The included and excluded covariate names are reported in
`diagnostics$spatial_iv_covariates` and
`diagnostics$excluded_spatial_iv_covariates`.

## Returned values

`PSPSAR()` returns coefficient and spillover tables, residuals,
HAC bandwidth information, and diagnostics. `SAR()` returns coefficient
tables, residuals, and diagnostics. See `?PSPSAR` and `?SAR` after
installation for the complete interface and return value definitions.

## References

Chen, Haiqiang, Yingxing Li, Xiaojun Song, and Yishu Wang (2026). *Profiled
Semiparametric GMM for a State-Varying SAR Coefficient*. Working paper.

Lan, Zhiqiang, Xiangfu Luo, Yishu Wang, and Guoyao Wu (2026). *Air Quality
and Spatial Spillovers in Electric Vehicle Charging: A State-Varying Spatial
Autoregressive Model*. Working paper.

Run `citation("PSPSAR")` after installation to retrieve these references in R.

## License and contact

The package is distributed under the [MIT License](LICENSE.md). For questions,
contact Yishu Wang at <yiswang@szu.edu.cn>.
