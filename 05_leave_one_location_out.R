suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

# Complete leave-one-location-out influence analysis for the registered DLM.
# Requires 00_gbd_analysis_core.R. The reference cumulative estimate is supplied
# directly in memory rather than read from a prior result file.

run_leave_one_location_out <- function(outcome, exposure, sdi, units, reference_cumulative) {
  prepared <- prepare_gbd_dlm_data(outcome, exposure, sdi, units)
  d <- prepared$d
  interaction_names <- prepared$interaction_names
  time_names <- prepared$time_names
  fixed_terms <- c("cause_ac", prepared$lag_names, interaction_names, "sdi_conditional", time_names)
  formula <- as.formula(paste("y_log ~", paste(fixed_terms, collapse = " + "), "+ (1 + time_std | location_factor)"))

  fit_one <- function(train) {
    fitted <- fit_lmer_with_fallback(formula, train)
    if (is.null(fitted$fit)) {
      return(list(status = "failed", cumulative = NA_real_, curve = rep(NA_real_, length(prepared$lag_values)), singular = NA, message = fitted$message))
    }
    fit <- fitted$fit
    curve_error <- ""
    curve <- tryCatch(
      extract_dlm_curve(fit, interaction_names, prepared$lag_transform, prepared$lag_spline),
      error = function(e) { curve_error <<- conditionMessage(e); rep(NA_real_, length(prepared$lag_values)) }
    )
    list(
      status = if (all(is.finite(curve))) "success" else "failed",
      cumulative = if (all(is.finite(curve))) mean(curve) else NA_real_,
      curve = curve,
      singular = isSingular(fit, tol = 1e-4),
      message = curve_error
    )
  }

  full <- fit_one(d)
  if (full$status != "success") stop("Full reference refit failed")
  location_ids <- sort(unique(d$location_id))
  if (length(location_ids) != 204L) stop("Expected 204 locations")
  unit_table <- assert_primary_units(units)

  records <- lapply(location_ids, function(location_id) {
    drop_id <- location_id
    result <- fit_one(d[location_id != drop_id])
    difference <- result$cumulative - reference_cumulative
    data.table(
      location_id = drop_id,
      location_name = if ("location_name" %in% names(unit_table)) as.character(unit_table$location_name[match(drop_id, unit_table$location_id)]) else NA_character_,
      removed_rows = sum(d$location_id == drop_id),
      fit_status = result$status,
      singular = result$singular,
      cumulative_delta = result$cumulative,
      difference_vs_reference = difference,
      absolute_difference_vs_reference = abs(difference),
      relative_difference_vs_reference = difference / reference_cumulative,
      message = result$message
    )
  })
  influence <- rbindlist(records, fill = TRUE)
  setorder(influence, -absolute_difference_vs_reference, location_id)

  list(
    full_refit_cumulative = full$cumulative,
    influence = influence,
    diagnostics = list(
      deletion_fit_count = nrow(influence),
      successful_deletion_fit_count = sum(influence$fit_status == "success"),
      singular_fit_count = sum(influence$singular %in% TRUE, na.rm = TRUE),
      all_cumulative_values_positive = all(influence[fit_status == "success", cumulative_delta] > 0),
      max_absolute_difference = max(influence[fit_status == "success", absolute_difference_vs_reference], na.rm = TRUE)
    )
  )
}
