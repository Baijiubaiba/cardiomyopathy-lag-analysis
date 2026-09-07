suppressPackageStartupMessages(library(data.table))

# Same-source temporal and reference-region blocked prediction stress tests.
# Requires 00_gbd_analysis_core.R to be loaded. The primary cumulative estimate
# is passed directly from the primary analysis result rather than read from disk.

run_blocked_validation <- function(outcome, exposure, sdi, units, source_dlm_cumulative) {
  prepared <- prepare_gbd_dlm_data(outcome, exposure, sdi, units)
  structures <- list(
    list(id = "random_intercept_slope", random_effect_term = "(1 + time_std | location_factor)")
  )
  result <- run_blocked_model_comparison(
    prepared = prepared,
    source_dlm_cumulative = source_dlm_cumulative,
    structures = structures,
    expected_regions = 21L
  )
  aggregate <- result$validation[, .(
    mean_rmse_log = mean(rmse_log),
    mean_mae_log = mean(mae_log),
    mean_bias_log = mean(bias_log)
  ), by = .(validation_type, model)]
  setorder(aggregate, validation_type, model)
  result$aggregate <- aggregate
  result
}
