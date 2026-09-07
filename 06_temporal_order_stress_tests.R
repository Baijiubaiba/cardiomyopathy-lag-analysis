suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(lme4)
})

# Temporal-order stress tests: past-aligned 0-3-year mean, future 1-4-year lead,
# and within-location random-time permutation. All data are passed in memory.

run_temporal_order_stress_tests <- function(
  outcome,
  exposure,
  sdi,
  units,
  bootstrap_replicates = 300L,
  bootstrap_seed = 20260904L,
  minimum_successful_bootstrap = 270L,
  permutation_replicates = 20L,
  permutation_seed = 20260903L,
  minimum_successful_permutations = 18L
) {
  units <- assert_primary_units(units)
  ids <- units$location_id
  analysis_years <- 2000:2019
  support_years <- 1996:2023
  source_years <- 1990:2023

  outcome <- as.data.table(copy(outcome))
  exposure <- as.data.table(copy(exposure))
  sdi <- as.data.table(copy(sdi))
  outcome[, c("location_id", "year", "cause_id", "measure_id", "metric_id", "age_id", "sex_id") := lapply(.SD, as.integer), .SDcols = c("location_id", "year", "cause_id", "measure_id", "metric_id", "age_id", "sex_id")]
  exposure[, c("location_id", "year", "rei_id", "measure_id", "metric_id", "age_id", "sex_id") := lapply(.SD, as.integer), .SDcols = c("location_id", "year", "rei_id", "measure_id", "metric_id", "age_id", "sex_id")]
  sdi[, c("location_id", "year_id", "age_group_id", "sex_id") := lapply(.SD, as.integer), .SDcols = c("location_id", "year_id", "age_group_id", "sex_id")]

  outcome <- outcome[cause_id %in% c(938L, 944L) & measure_id == 2L & metric_id == 3L & age_id == 22L & sex_id == 3L & year %in% source_years & location_id %in% ids]
  exposure <- exposure[rei_id == 102L & measure_id == 29L & metric_id == 3L & age_id == 22L & sex_id == 3L & year %in% source_years & location_id %in% ids]
  sdi <- sdi[covariate_name_short == "sdi" & age_group_id == 22L & sex_id == 3L & year_id %in% source_years & location_id %in% ids]
  if (nrow(outcome) != 204L * 34L * 2L || nrow(exposure) != 204L * 34L || nrow(sdi) != 204L * 34L) stop("Unexpected temporal-stress-test source coverage")

  outcome <- outcome[year %in% analysis_years]
  exposure_lookup <- exposure[, .(location_id, year, exposure = as.numeric(val))]
  sdi_lookup <- sdi[, .(location_id, year = year_id, sdi = as.numeric(mean_value))]
  exposure_stats <- exposure_lookup[year %in% support_years, .(mean = mean(exposure), sd = sd(exposure))]
  sdi_stats <- sdi_lookup[year %in% analysis_years, .(mean = mean(sdi), sd = sd(sdi))]
  if (exposure_stats$sd <= 0 || sdi_stats$sd <= 0) stop("Standardization scales must be positive")

  base_keys <- unique(outcome[, .(location_id, year)])
  build_feature <- function(lookup, direction, offsets) {
    z <- copy(base_keys)
    value_names <- character()
    for (offset in offsets) {
      if (direction == "past") {
        tmp <- lookup[, .(location_id, outcome_year = year + offset, x = exposure)]
      } else if (direction == "future") {
        tmp <- lookup[, .(location_id, outcome_year = year - offset, x = exposure)]
      } else stop("Unknown temporal direction")
      value_name <- paste0("v", offset)
      setnames(tmp, "x", value_name)
      z <- merge(z, tmp, by.x = c("location_id", "year"), by.y = c("location_id", "outcome_year"), all.x = TRUE, sort = FALSE)
      value_names <- c(value_names, value_name)
    }
    if (anyNA(z[, ..value_names])) stop("Incomplete temporal exposure feature")
    raw_mean <- rowMeans(z[, ..value_names])
    data.table(location_id = z$location_id, year = z$year, feature = (raw_mean - exposure_stats$mean) / exposure_stats$sd)
  }

  past <- build_feature(exposure_lookup, "past", 0:3)
  future <- build_feature(exposure_lookup, "future", 1:4)
  setnames(past, "feature", "past_mean_0_3")
  setnames(future, "feature", "future_lead_mean_1_4")
  outcome <- merge(outcome, past, by = c("location_id", "year"), all.x = TRUE, sort = FALSE)
  outcome <- merge(outcome, future, by = c("location_id", "year"), all.x = TRUE, sort = FALSE)
  outcome <- merge(outcome, sdi_lookup, by = c("location_id", "year"), all.x = TRUE, sort = FALSE)
  if (nrow(outcome) != 204L * length(analysis_years) * 2L || anyNA(outcome[, .(past_mean_0_3, future_lead_mean_1_4, sdi)])) stop("Incomplete temporal model table")

  outcome[, time_idx := match(year, analysis_years)]
  outcome[, time_std := 2 * (time_idx - 1) / (length(analysis_years) - 1) - 1]
  outcome[, cause_ac := as.integer(cause_id == 938L)]
  outcome[, y_log := log(as.numeric(val))]
  outcome[, sdi_std := (sdi - sdi_stats$mean) / sdi_stats$sd]
  calendar_spline <- ns(analysis_years, df = 3, Boundary.knots = c(2000, 2019))
  time_features <- calendar_spline[outcome$time_idx, , drop = FALSE]
  sdi_conditional <- as.numeric(residuals(lm(outcome$sdi_std ~ time_features)))

  d <- data.table(
    location_id = outcome$location_id,
    year = outcome$year,
    y_log = outcome$y_log,
    cause_ac = outcome$cause_ac,
    sdi_conditional = sdi_conditional,
    time_std = outcome$time_std,
    past_mean_0_3 = outcome$past_mean_0_3,
    future_lead_mean_1_4 = outcome$future_lead_mean_1_4
  )
  d[, location_factor := factor(as.character(location_id), levels = as.character(sort(unique(location_id))))]
  time_names <- paste0("time_q", seq_len(ncol(time_features)))
  for (j in seq_len(ncol(time_features))) d[[time_names[j]]] <- time_features[, j]

  fit_one <- function(data, feature_name, group_name) {
    interaction_name <- paste0("cause_", feature_name)
    data[[interaction_name]] <- data$cause_ac * data[[feature_name]]
    fixed_terms <- c("cause_ac", feature_name, interaction_name, "sdi_conditional", time_names)
    formula <- as.formula(paste("y_log ~", paste(fixed_terms, collapse = " + "), "+ (1 + time_std |", group_name, ")"))
    fitted <- fit_lmer_with_fallback(formula, data)
    if (is.null(fitted$fit)) return(list(status = "failed", estimate = NA_real_, ratio = NA_real_, singular = NA, message = fitted$message))
    fit <- fitted$fit
    fixed <- fixef(fit)
    if (!(interaction_name %in% names(fixed)) || !is.finite(fixed[[interaction_name]])) return(list(status = "failed", estimate = NA_real_, ratio = NA_real_, singular = isSingular(fit, tol = 1e-4), message = "Missing temporal interaction coefficient"))
    estimate <- as.numeric(fixed[[interaction_name]])
    list(status = "success", estimate = estimate, ratio = exp(estimate), singular = isSingular(fit, tol = 1e-4), message = "")
  }

  temporal_ids <- c("past_mean_0_3", "future_lead_mean_1_4")
  full_results <- setNames(lapply(temporal_ids, function(x) fit_one(copy(d), x, "location_factor")), temporal_ids)
  full <- rbindlist(lapply(temporal_ids, function(x) {
    r <- full_results[[x]]
    data.table(control = x, fit_status = r$status, singular = r$singular, estimate_log_rate_contrast = r$estimate, model_implied_ratio = r$ratio, message = r$message)
  }))

  unit_rows <- split(d, d$location_id)
  unit_ids <- names(unit_rows)
  set.seed(bootstrap_seed)
  bootstrap_records <- vector("list", bootstrap_replicates * length(temporal_ids))
  k <- 0L
  for (b in seq_len(bootstrap_replicates)) {
    sampled <- sample(unit_ids, length(unit_ids), replace = TRUE)
    boot <- rbindlist(lapply(seq_along(sampled), function(j) {
      z <- copy(unit_rows[[sampled[j]]])
      z[, boot_loc := paste0(sampled[j], "_", j)]
      z
    }), fill = TRUE)
    boot[, boot_loc := factor(boot_loc)]
    for (feature_name in temporal_ids) {
      r <- fit_one(copy(boot), feature_name, "boot_loc")
      k <- k + 1L
      bootstrap_records[[k]] <- data.table(replicate = b, control = feature_name, status = r$status, estimate = r$estimate, ratio = r$ratio, singular = r$singular, message = r$message)
    }
  }
  bootstrap <- rbindlist(bootstrap_records, fill = TRUE)
  success <- bootstrap[status == "success" & is.finite(estimate), .N, by = control]
  if (nrow(success[N >= minimum_successful_bootstrap]) != length(temporal_ids)) stop("Temporal bootstrap success below registered minimum")

  bootstrap_summary <- bootstrap[status == "success" & is.finite(estimate), .(
    bootstrap_success = .N,
    bootstrap_lower = quantile(estimate, 0.025, names = FALSE),
    bootstrap_median = quantile(estimate, 0.5, names = FALSE),
    bootstrap_upper = quantile(estimate, 0.975, names = FALSE),
    bootstrap_singular = sum(singular %in% TRUE, na.rm = TRUE)
  ), by = control]

  paired <- dcast(bootstrap[status == "success" & is.finite(estimate)], replicate ~ control, value.var = "estimate")
  paired <- paired[is.finite(past_mean_0_3) & is.finite(future_lead_mean_1_4)]
  if (nrow(paired) < minimum_successful_bootstrap) stop("Paired temporal bootstrap success below registered minimum")
  difference <- paired$future_lead_mean_1_4 - paired$past_mean_0_3
  comparison <- data.table(
    comparison = "future_lead_minus_past_aligned",
    full_difference_log_rate_contrast = full_results$future_lead_mean_1_4$estimate - full_results$past_mean_0_3$estimate,
    paired_bootstrap_success = length(difference),
    paired_bootstrap_lower = quantile(difference, 0.025, names = FALSE),
    paired_bootstrap_median = quantile(difference, 0.5, names = FALSE),
    paired_bootstrap_upper = quantile(difference, 0.975, names = FALSE)
  )

  exposure_support <- exposure_lookup[year %in% support_years]
  permutation_records <- vector("list", permutation_replicates)
  for (r in seq_len(permutation_replicates)) {
    set.seed(permutation_seed + r - 1L)
    perm_lookup <- copy(exposure_support)
    perm_lookup[, exposure := sample(exposure), by = location_id]
    feature <- build_feature(perm_lookup, "past", 0:3)
    setnames(feature, "feature", "random_time_permutation_mean_0_3")
    perm <- merge(copy(d), feature, by = c("location_id", "year"), all.x = TRUE, sort = FALSE)
    perm[, location_factor := factor(as.character(location_id), levels = as.character(sort(unique(location_id))))]
    fit <- fit_one(perm, "random_time_permutation_mean_0_3", "location_factor")
    permutation_records[[r]] <- data.table(replicate = r, fit_status = fit$status, singular = fit$singular, estimate_log_rate_contrast = fit$estimate, model_implied_ratio = fit$ratio, message = fit$message)
  }
  permutation <- rbindlist(permutation_records, fill = TRUE)
  valid_permutation <- permutation[fit_status == "success" & is.finite(estimate_log_rate_contrast)]
  if (nrow(valid_permutation) < minimum_successful_permutations) stop("Permutation success below registered minimum")
  permutation_summary <- data.table(
    control = "random_time_permutation_mean_0_3",
    permutation_success = nrow(valid_permutation),
    estimate_mean = mean(valid_permutation$estimate_log_rate_contrast),
    empirical_lower = quantile(valid_permutation$estimate_log_rate_contrast, 0.025, names = FALSE),
    empirical_median = quantile(valid_permutation$estimate_log_rate_contrast, 0.5, names = FALSE),
    empirical_upper = quantile(valid_permutation$estimate_log_rate_contrast, 0.975, names = FALSE),
    singular_count = sum(valid_permutation$singular %in% TRUE, na.rm = TRUE)
  )

  list(
    full = full,
    bootstrap = bootstrap,
    bootstrap_summary = bootstrap_summary,
    comparison = comparison,
    permutation = permutation,
    permutation_summary = permutation_summary
  )
}
