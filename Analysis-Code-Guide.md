# Analysis code use guide

The accompanying Clean-Analysis-Code.zip contains eleven analysis source files, ten in R and one in Python. This guide is outside that archive. Data are acquired separately from the providers identified in the manuscript. The functions receive data frames or prepared in-memory objects and return analysis objects; they do not read private project paths or write result files.

## Input definitions

The GBD outcome frame requires location_id, year, cause_id, measure_id, metric_id, age_id, sex_id, val, lower and upper. The retained outcome codes are AC 938 and OC 944; DALYs measure 2, Rate metric 3, All ages 22 and Both sexes 3. Input support spans 1990–2023; the main outcome window is 2000–2023.

The GBD alcohol exposure frame requires location_id, year, rei_id, measure_id, metric_id, age_id, sex_id and val. The retained exposure is High alcohol use SEV, rei_id 102, measure 29, metric 3, age 22 and sex 3.

The SDI frame requires location_id, year_id, covariate_name_short, age_group_id, sex_id and mean_value. The selected covariate is sdi, age 22 and sex 3. The units frame identifies the 204 non-overlapping countries/territories, using location_id and optionally primary_eligible. Region-block and region-balanced analyses additionally require reference_parent_id_2021, the compatibility grouping documented in Supplementary Table S1. Unit names and coverage are reported in the supplementary tables.

The WHO analysis additionally requires a country_crosswalk with gbd_location_id, gbd_location_name, iso3 and who_code_present. Its WHO frame requires IND_CODE, DIM_TIME, DIM_TIME_TYPE, DIM_GEO_CODE_M49, DIM_GEO_CODE_TYPE, DIM_PUBLISH_STATE_CODE, GEO_NAME_SHORT, DIM_SEX and RATE_PER_CAPITA_N. The code filters SA_0000001688, published country-year records and TOTAL sex, retaining 188 complete 2000–2022 trajectories.

## Source-file order and entry points

1. Load 00_gbd_analysis_core.R first. prepare_gbd_dlm_data(outcome, exposure, sdi, units) constructs the original standardized lag features and returns the prepared object. It checks coverage and pairing prerequisites.

2. Load 01_primary_cumulative_model.R for run_primary_cumulative_model(). This returns the original Bayesian posterior summaries and country-bootstrap objects. The returned delta and summary objects supply the primary curve and cumulative estimates needed by comparison functions.

3. Load 02_blocked_validation.R and 03_random_effects_sensitivity.R for the corresponding run_blocked_validation() and run_random_effects_sensitivity() functions. Both receive outcome, exposure, sdi, units and source_dlm_cumulative.

4. The Python file 04_effect_scale.py exposes compute_effect_scale(). It receives the exposure frame, primary location IDs and the primary/bootstrap cumulative summaries, returning model-implied scale translations. It is independent of the R loading sequence.

5. Load 05_leave_one_location_out.R for run_leave_one_location_out(), passing the original full-data maximum-likelihood reference cumulative estimate. Load 06_temporal_order_stress_tests.R and 07_weighting_sensitivity.R for the original temporal-order and weighting checks. Their seeds and replication defaults are specified in the function signatures.

6. Load 08_lag_shape_sensitivity.R for run_lag_shape_sensitivity(), supplying the original primary curve and cumulative estimate. Load 09_who_exposure_triangulation.R for run_who_exposure_triangulation(), with the additional WHO frame and crosswalk defined above.

7. Load 10_differential_within_between.R for analyse_differential_within_between(prepared). The prepared object comes from step 1. Its defaults retain 300 matched country resamples, seed 20260906 and a minimum of 285 successful replicates. It returns the four differential-adjustment summaries, the between/within cumulative summaries and the paired within-minus-between comparison.

## Dependencies and verification scope

The R sources use data.table, splines and lme4; the original Bayesian module also requires JAGS, rjags and coda. The Python scale translation uses NumPy and pandas. The targeted secondary analysis was executed using R 4.6.0, data.table 1.18.4 and lme4 2.0.1. These version statements describe that recorded run, not a reconstructed historical environment for every original analysis.

The secondary source implementation was numerically compared with the recorded analysis summaries and matched country-bootstrap intervals. This is computational verification of the specified implementation, not independent clinical replication or a demonstrated third-party end-to-end rerun of all historical modules. The normalized cumulative contrast is the average of eleven lag coefficients. Between and within estimates are different estimands, not contribution shares of the pooled coefficient.
