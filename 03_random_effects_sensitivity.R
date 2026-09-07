suppressPackageStartupMessages(library(data.table))

# Registered random-effects structure sensitivity. Requires the shared core.
# Both random-intercept-plus-time-slope and random-intercept-only structures are
# fitted under the same four exposure representations and blocked-validation design.

run_random_effects_sensitivity <- function(outcome, exposure, sdi, units, source_dlm_cumulative) {
  prepared <- prepare_gbd_dlm_data(outcome, exposure, sdi, units)
  structures <- list(
    list(id = "random_intercept_slope", random_effect_term = "(1 + time_std | location_factor)"),
    list(id = "random_intercept_only", random_effect_term = "(1 | location_factor)")
  )
  result <- run_blocked_model_comparison(
    prepared = prepared,
    source_dlm_cumulative = source_dlm_cumulative,
    structures = structures,
    expected_regions = 21L
  )
  result$dlm_only <- result$validation[model == "dlm_df4"]
  result$dlm_structure_summary <- result$dlm_only[, .(
    blocked_fit_count = .N,
    successful_fit_count = sum(fit_status == "success"),
    singular_fit_count = sum(singular %in% TRUE, na.rm = TRUE),
    mean_rmse_log = mean(rmse_log),
    all_cumulative_positive = all(cumulative_delta > 0)
  ), by = .(structure, validation_type)]
  result
}
