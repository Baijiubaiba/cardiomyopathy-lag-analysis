# Cardiomyopathy lag-association analysis

Analysis source for alcohol-related and other cardiomyopathy burden using GBD 2023 population estimates, with distributed-lag models and differential within/between-location analyses.

## Scope

The analysis covers 204 non-overlapping countries/territories. The outcome window is 2000-2023, with 1990-2023 exposure support for lags 0-10. Outcomes are all-age, both-sex DALY rates for alcoholic cardiomyopathy (AC) and other cardiomyopathy (OC).

This repository contains eleven analysis modules (ten R files and one Python file). Functions accept data frames or prepared objects in memory and return analysis objects. Data, manuscript files, generated results and internal project records are not distributed.

## Use

See [Analysis-Code-Guide.md](Analysis-Code-Guide.md) for required columns, filtering codes, dependencies, function arguments and execution order. Obtain source data separately under the respective providers' access terms. Exact country selection, regional grouping and WHO crosswalk must be supplied as described in the guide; these reference tables are not bundled.

Run from this directory after constructing the four input data frames:

```r
source("00_gbd_analysis_core.R")
prepared <- prepare_gbd_dlm_data(outcome, exposure, sdi, units)

# Original hierarchical distributed-lag analysis (requires JAGS).
source("01_primary_cumulative_model.R")
primary <- run_primary_cumulative_model(outcome, exposure, sdi, units)

# Subsequent differential-adjustment and within/between analysis.
source("10_differential_within_between.R")
secondary <- analyse_differential_within_between(prepared)
```

The original analysis precedes the subsequent differential analysis. The latter retains the original preprocessing and does not refit the Bayesian model. Other modules provide blocked validation, structural sensitivity, effect-scale translation, location deletion, temporal-order checks, weighting, lag-shape sensitivity and WHO exposure-source triangulation.

## Dependencies

R: `data.table`, `splines` (included with R), `lme4`; the Bayesian module also uses JAGS, `rjags` and `coda`. Python: NumPy and pandas. Package-version information for the recorded secondary run is in the guide; it is not a complete environment lock for all historical analyses.

## Interpretation and reproducibility

The normalized cumulative contrast is the average of eleven lag coefficients. Within-location and between-location estimates are distinct estimands, not shares of the pooled coefficient. Estimates describe associations among population-level estimates; they do not establish causality, clinical risk, biological lag windows or independent outcome validation. WHO comparisons change the exposure source while retaining GBD outcomes.

This is a source-code release, not a self-contained data-and-environment reproduction bundle. Parsing and release-content checks do not certify an end-to-end rerun. Historical numerical checks of the secondary implementation are described in the guide.

## Release

Initial source release: 2026-09-08. The eleven analysis modules are unchanged from the selected submission source package. This publication does not constitute journal submission or final manuscript approval.

No software license has been selected for this release. Public availability alone does not grant an open-source license. Data-provider rights and terms remain separate.
