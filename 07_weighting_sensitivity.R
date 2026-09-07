suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

# Equal-unit versus region-balanced weighting sensitivity for the registered DLM.
# Population weighting is intentionally not computed because the registered study
# did not contain the required traceable GBD 2023 country-year population input.

run_weighting_sensitivity <- function(
  outcome,
  exposure,
  sdi,
  units,
  bootstrap_replicates = 300L,
  seed = 20260905L,
  minimum_successful_bootstrap = 285L
) {
  prepared <- prepare_gbd_dlm_data(outcome, exposure, sdi, units)
  d <- copy(prepared$d)
  unit_table <- prepared$units
  if (!"reference_parent_id_2021" %in% names(unit_table)) stop("units must contain reference_parent_id_2021")

  location_region <- data.table(
    location_id = unit_table$location_id,
    region_id = as.character(unit_table$reference_parent_id_2021)
  )
  if (uniqueN(location_region$region_id) != 21L || anyNA(location_region$region_id)) stop("Expected 21 reference regions")
  region_counts <- location_region[, .(n_locations = .N), by = region_id]
  n_regions <- uniqueN(location_region$region_id)
  location_region <- merge(location_region, region_counts, by = "region_id", all.x = TRUE, sort = FALSE)
  location_region[, region_balanced_weight := nrow(location_region) / (n_regions * n_locations)]
  d <- merge(d, location_region[, .(location_id, region_id, region_balanced_weight)], by = "location_id", all.x = TRUE, sort = FALSE)
  d[, equal_unit_weight := 1]
  if (anyNA(d$region_balanced_weight) || any(d$region_balanced_weight <= 0)) stop("Invalid region-balanced weights")

  fixed_terms <- c("cause_ac", prepared$lag_names, prepared$interaction_names, "sdi_conditional", prepared$time_names)
  formula_for_group <- function(group_name) {
    as.formula(paste("y_log ~", paste(fixed_terms, collapse = " + "), "+ (1 + time_std |", group_name, ")"))
  }
  fit_one <- function(data, weight_name, group_name) {
    fit_data <- copy(data)
    weights <- as.numeric(fit_data[[weight_name]])
    if (any(!is.finite(weights)) || any(weights <= 0)) stop("Invalid analysis weights")
    fitted <- fit_lmer_with_fallback(formula_for_group(group_name), fit_data, weights = weights)
    if (is.null(fitted$fit)) return(list(status = "failed", cumulative = NA_real_, curve = rep(NA_real_, length(prepared$lag_values)), singular = NA, message = fitted$message))
    fit <- fitted$fit
    curve_error <- ""
    curve <- tryCatch(
      extract_dlm_curve(fit, prepared$interaction_names, prepared$lag_transform, prepared$lag_spline),
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

  schemes <- c("equal_unit_weight", "region_balanced_weight")
  full <- rbindlist(lapply(schemes, function(scheme) {
    result <- fit_one(d, scheme, "location_factor")
    data.table(
      scheme = scheme,
      fit_status = result$status,
      singular = result$singular,
      cumulative_delta = result$cumulative,
      message = result$message
    )
  }))

  unit_rows <- split(d, d$location_id)
  unit_ids <- names(unit_rows)
  set.seed(seed)
  records <- vector("list", bootstrap_replicates * length(schemes))
  k <- 0L
  for (b in seq_len(bootstrap_replicates)) {
    sampled <- sample(unit_ids, length(unit_ids), replace = TRUE)
    boot <- rbindlist(lapply(seq_along(sampled), function(j) {
      z <- copy(unit_rows[[sampled[j]]])
      z[, boot_loc := paste0(sampled[j], "_", j)]
      z
    }), fill = TRUE)
    boot[, boot_loc := factor(boot_loc)]
    for (scheme in schemes) {
      result <- fit_one(boot, scheme, "boot_loc")
      k <- k + 1L
      records[[k]] <- data.table(
        replicate = b,
        scheme = scheme,
        status = result$status,
        cumulative_delta = result$cumulative,
        singular = result$singular,
        message = result$message
      )
    }
  }
  bootstrap <- rbindlist(records, fill = TRUE)
  summary <- bootstrap[status == "success" & is.finite(cumulative_delta), .(
    bootstrap_success = .N,
    bootstrap_lower = quantile(cumulative_delta, 0.025, names = FALSE),
    bootstrap_median = quantile(cumulative_delta, 0.5, names = FALSE),
    bootstrap_upper = quantile(cumulative_delta, 0.975, names = FALSE),
    bootstrap_singular = sum(singular %in% TRUE, na.rm = TRUE)
  ), by = scheme]
  if (nrow(summary[bootstrap_success >= minimum_successful_bootstrap]) != length(schemes)) stop("Weighting bootstrap success below registered minimum")
  summary <- merge(full, summary, by = "scheme", all.x = TRUE, sort = FALSE)
  summary <- rbind(
    summary,
    data.table(
      scheme = "population_weight",
      fit_status = "not_assessable",
      singular = NA,
      cumulative_delta = NA_real_,
      message = "Traceable GBD 2023 country-year population input was not available in the registered study"
    ),
    fill = TRUE
  )
  list(summary = summary, bootstrap = bootstrap, location_weights = location_region)
}
