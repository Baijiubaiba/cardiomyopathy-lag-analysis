suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(lme4)
})

# Shared in-memory analysis helpers for the registered GBD 2023 analyses.
# No file-system, project-state, reporting, figure, or table-building operations
# are performed here. Input data frames are supplied by the analyst.

assert_primary_units <- function(units, expected_n = 204L) {
  units <- as.data.table(copy(units))
  if (!"location_id" %in% names(units)) stop("units must contain location_id")
  units[, location_id := as.integer(location_id)]
  if ("primary_eligible" %in% names(units)) units <- units[primary_eligible %in% TRUE]
  if (nrow(units) != expected_n || uniqueN(units$location_id) != expected_n) {
    stop(sprintf("Expected %d unique primary units", expected_n))
  }
  units
}

prepare_gbd_dlm_data <- function(
  outcome,
  exposure,
  sdi,
  units,
  analysis_years = 2000:2023,
  source_years = 1990:2023,
  lag_values = 0:10,
  exposure_standardization_years = analysis_years
) {
  units <- assert_primary_units(units)
  primary_ids <- units$location_id
  outcome <- as.data.table(copy(outcome))
  exposure <- as.data.table(copy(exposure))
  sdi <- as.data.table(copy(sdi))

  required_outcome <- c("location_id", "year", "cause_id", "measure_id", "metric_id", "age_id", "sex_id", "val", "lower", "upper")
  required_exposure <- c("location_id", "year", "rei_id", "measure_id", "metric_id", "age_id", "sex_id", "val")
  required_sdi <- c("location_id", "year_id", "covariate_name_short", "age_group_id", "sex_id", "mean_value")
  if (!all(required_outcome %in% names(outcome))) stop("outcome is missing required GBD fields")
  if (!all(required_exposure %in% names(exposure))) stop("exposure is missing required GBD fields")
  if (!all(required_sdi %in% names(sdi))) stop("sdi is missing required GBD fields")

  outcome[, c("location_id", "year", "cause_id", "measure_id", "metric_id", "age_id", "sex_id") :=
    lapply(.SD, as.integer), .SDcols = c("location_id", "year", "cause_id", "measure_id", "metric_id", "age_id", "sex_id")]
  exposure[, c("location_id", "year", "rei_id", "measure_id", "metric_id", "age_id", "sex_id") :=
    lapply(.SD, as.integer), .SDcols = c("location_id", "year", "rei_id", "measure_id", "metric_id", "age_id", "sex_id")]
  sdi[, c("location_id", "year_id", "age_group_id", "sex_id") :=
    lapply(.SD, as.integer), .SDcols = c("location_id", "year_id", "age_group_id", "sex_id")]

  outcome <- outcome[
    cause_id %in% c(938L, 944L) & measure_id == 2L & metric_id == 3L &
      age_id == 22L & sex_id == 3L & year %in% source_years & location_id %in% primary_ids
  ]
  exposure <- exposure[
    rei_id == 102L & measure_id == 29L & metric_id == 3L & age_id == 22L &
      sex_id == 3L & year %in% source_years & location_id %in% primary_ids
  ]
  sdi <- sdi[
    covariate_name_short == "sdi" & age_group_id == 22L & sex_id == 3L &
      year_id %in% source_years & location_id %in% primary_ids
  ]

  expected_source_years <- length(source_years)
  if (nrow(outcome) != 204L * expected_source_years * 2L) stop("Unexpected AC/OC source coverage")
  if (nrow(exposure) != 204L * expected_source_years) stop("Unexpected alcohol-exposure source coverage")
  if (nrow(sdi) != 204L * expected_source_years) stop("Unexpected SDI source coverage")
  if (outcome[, .N, by = .(location_id, year, cause_id)][, any(N != 1L)]) stop("Duplicate outcome keys")
  if (exposure[, .N, by = .(location_id, year, rei_id)][, any(N != 1L)]) stop("Duplicate exposure keys")
  if (sdi[, .N, by = .(location_id, year_id)][, any(N != 1L)]) stop("Duplicate SDI keys")
  if (any(!is.finite(outcome$val)) || any(outcome$val <= 0)) stop("Outcome point estimates must be positive and finite")

  outcome <- outcome[year %in% analysis_years]
  exposure_lookup <- exposure[, .(location_id, year, exposure = as.numeric(val))]
  sdi_lookup <- sdi[, .(location_id, year = year_id, sdi = as.numeric(mean_value))]
  exposure_stats <- exposure_lookup[year %in% exposure_standardization_years, .(mean = mean(exposure), sd = sd(exposure))]
  sdi_stats <- sdi_lookup[year %in% analysis_years, .(mean = mean(sdi), sd = sd(sdi))]
  if (exposure_stats$sd <= 0 || sdi_stats$sd <= 0) stop("Standardization scales must be positive")

  base <- unique(outcome[, .(location_id, year)])
  lag_cols <- paste0("x", lag_values)
  for (lag in lag_values) {
    tmp <- exposure_lookup[, .(location_id, outcome_year = year + lag, x = exposure)]
    setnames(tmp, "x", paste0("x", lag))
    base <- merge(base, tmp, by.x = c("location_id", "year"), by.y = c("location_id", "outcome_year"), all.x = TRUE, sort = FALSE)
  }
  base <- merge(base, sdi_lookup, by = c("location_id", "year"), all.x = TRUE, sort = FALSE)
  if (anyNA(base[, ..lag_cols]) || anyNA(base$sdi)) stop("Missing lagged exposure or SDI")
  base[, (lag_cols) := lapply(.SD, function(z) (z - exposure_stats$mean) / exposure_stats$sd), .SDcols = lag_cols]
  base[, sdi_std := (sdi - sdi_stats$mean) / sdi_stats$sd]

  outcome <- merge(outcome, base[, c("location_id", "year", lag_cols, "sdi_std"), with = FALSE], by = c("location_id", "year"), all.x = TRUE, sort = FALSE)
  if (nrow(outcome) != 204L * length(analysis_years) * 2L) stop("Unexpected analysis-table row count")
  if (anyNA(outcome[, c(lag_cols, "sdi_std"), with = FALSE])) stop("Incomplete analysis table")

  outcome[, time_idx := match(year, analysis_years)]
  outcome[, time_std := 2 * (time_idx - 1) / (length(analysis_years) - 1) - 1]
  outcome[, cause_ac := as.integer(cause_id == 938L)]
  outcome[, y_log := log(as.numeric(val))]
  outcome[, location_factor := factor(as.character(location_id), levels = as.character(sort(unique(location_id))))]

  lag_spline <- ns(lag_values, df = 4, knots = c(3, 7), Boundary.knots = c(0, 10), intercept = TRUE)
  lag_input <- as.matrix(outcome[, ..lag_cols])
  lag_features_raw <- lag_input %*% lag_spline
  lag_center <- colMeans(lag_features_raw)
  lag_qr <- qr(sweep(lag_features_raw, 2L, lag_center, "-"), LAPACK = FALSE)
  lag_scale <- sqrt(nrow(lag_features_raw) - 1)
  lag_features <- qr.Q(lag_qr) * lag_scale
  lag_transform <- qr.R(lag_qr) / lag_scale

  calendar_spline <- ns(analysis_years, df = 3, Boundary.knots = c(min(analysis_years), max(analysis_years)))
  time_features <- calendar_spline[outcome$time_idx, , drop = FALSE]
  sdi_conditional <- as.numeric(residuals(lm(outcome$sdi_std ~ time_features)))

  d <- data.table(
    location_id = outcome$location_id,
    year = outcome$year,
    y_log = outcome$y_log,
    cause_ac = outcome$cause_ac,
    sdi_conditional = sdi_conditional,
    time_std = outcome$time_std,
    location_factor = outcome$location_factor
  )
  if ("reference_parent_id_2021" %in% names(units)) {
    d[, region_id := as.character(units$reference_parent_id_2021[match(location_id, units$location_id)])]
  }
  for (w in c(0L, 5L, 10L)) {
    mean_name <- paste0("exposure_mean_", w)
    d[[mean_name]] <- rowMeans(outcome[, paste0("x", 0:w), with = FALSE])
    d[[paste0("cause_", mean_name)]] <- d$cause_ac * d[[mean_name]]
  }

  lag_names <- paste0("dlm_q", seq_len(ncol(lag_features)))
  interaction_names <- paste0("cause_dlm_q", seq_len(ncol(lag_features)))
  time_names <- paste0("time_q", seq_len(ncol(time_features)))
  for (k in seq_len(ncol(lag_features))) d[[lag_names[k]]] <- lag_features[, k]
  for (k in seq_len(ncol(lag_features))) d[[interaction_names[k]]] <- d$cause_ac * d[[lag_names[k]]]
  for (j in seq_len(ncol(time_features))) d[[time_names[j]]] <- time_features[, j]

  list(
    units = units,
    outcome = outcome,
    d = d,
    exposure_lookup = exposure_lookup,
    sdi_lookup = sdi_lookup,
    exposure_stats = exposure_stats,
    sdi_stats = sdi_stats,
    analysis_years = analysis_years,
    lag_values = lag_values,
    lag_cols = lag_cols,
    lag_spline = lag_spline,
    lag_transform = lag_transform,
    lag_features = lag_features,
    time_features = time_features,
    sdi_conditional = sdi_conditional,
    lag_names = lag_names,
    interaction_names = interaction_names,
    time_names = time_names
  )
}

fit_lmer_with_fallback <- function(formula, data, weights = NULL, maxfun = 100000L) {
  fit_once <- function(optimizer) {
    args <- list(
      formula = formula,
      data = data,
      REML = FALSE,
      control = lmerControl(
        optimizer = optimizer,
        optCtrl = list(maxfun = maxfun),
        check.rankX = "silent.drop.cols",
        check.conv.singular = "ignore"
      )
    )
    if (!is.null(weights)) args$weights <- weights
    suppressWarnings(do.call(lmer, args))
  }
  first_error <- ""
  fit <- tryCatch(fit_once("bobyqa"), error = function(e) { first_error <<- conditionMessage(e); NULL })
  if (is.null(fit)) {
    second_error <- ""
    fit <- tryCatch(fit_once("nloptwrap"), error = function(e) { second_error <<- conditionMessage(e); NULL })
    if (is.null(fit)) first_error <- paste0(first_error, " | fallback: ", second_error)
  }
  list(fit = fit, message = first_error)
}

extract_dlm_curve <- function(fit, interaction_names, lag_transform, lag_spline) {
  fixed <- fixef(fit)
  if (!all(interaction_names %in% names(fixed))) stop("Missing DLM interaction coefficients")
  gamma <- as.numeric(fixed[interaction_names])
  if (any(!is.finite(gamma))) stop("Non-finite DLM interaction coefficients")
  delta_original <- as.numeric(solve(lag_transform, gamma))
  as.numeric(lag_spline %*% delta_original)
}

dlm_model_definitions <- function() {
  list(
    list(id = "dlm_df4", label = "DLM df=4", kind = "dlm", window_end = NA_integer_, normalization_divisor = 11),
    list(id = "current_exposure", label = "Current exposure", kind = "mean", window_end = 0L, normalization_divisor = 1),
    list(id = "mean_5y", label = "5-year mean", kind = "mean", window_end = 5L, normalization_divisor = 6),
    list(id = "mean_10y", label = "10-year mean", kind = "mean", window_end = 10L, normalization_divisor = 11)
  )
}

run_blocked_model_comparison <- function(
  prepared,
  source_dlm_cumulative,
  structures,
  temporal_holdouts = list(
    list(id = "time_train_through_2015", train_end = 2015L, test_start = 2016L, test_end = 2019L),
    list(id = "time_train_through_2019", train_end = 2019L, test_start = 2020L, test_end = 2023L)
  ),
  expected_regions = 21L
) {
  d <- prepared$d
  if (!"region_id" %in% names(d)) stop("units must contain reference_parent_id_2021 for blocked-region analyses")
  if (!is.finite(source_dlm_cumulative)) stop("source_dlm_cumulative must be finite")
  models <- dlm_model_definitions()
  lag_values <- prepared$lag_values
  lag_names <- prepared$lag_names
  lag_interaction_names <- prepared$interaction_names
  time_names <- prepared$time_names

  model_features <- function(info) if (identical(info$kind, "dlm")) lag_names else paste0("exposure_mean_", info$window_end)
  model_interactions <- function(info) paste0("cause_", model_features(info))
  model_formula <- function(info, structure) {
    fixed_terms <- c("cause_ac", model_features(info), model_interactions(info), "sdi_conditional", time_names)
    as.formula(paste("y_log ~", paste(fixed_terms, collapse = " + "), "+", structure$random_effect_term))
  }
  extract_delta <- function(info, fit) {
    names_i <- model_interactions(info)
    fixed <- fixef(fit)
    if (!all(names_i %in% names(fixed))) stop("Missing model interaction coefficient")
    gamma <- as.numeric(fixed[names_i])
    if (identical(info$kind, "dlm")) {
      delta_original <- as.numeric(solve(prepared$lag_transform, gamma))
      as.numeric(prepared$lag_spline %*% delta_original)
    } else {
      rep(gamma[1] / info$normalization_divisor, length(lag_values))
    }
  }
  fit_one <- function(info, structure, train_idx, test_idx, include_random_effects) {
    fitted <- fit_lmer_with_fallback(model_formula(info, structure), d[train_idx])
    if (is.null(fitted$fit)) {
      return(list(status = "failed", prediction = rep(NA_real_, length(test_idx)), curve = rep(NA_real_, length(lag_values)), cumulative = NA_real_, singular = NA, message = fitted$message))
    }
    fit <- fitted$fit
    prediction_error <- ""
    prediction <- tryCatch(
      as.numeric(predict(fit, newdata = d[test_idx], re.form = if (include_random_effects) NULL else NA, allow.new.levels = TRUE)),
      error = function(e) { prediction_error <<- conditionMessage(e); rep(NA_real_, length(test_idx)) }
    )
    curve_error <- ""
    curve <- tryCatch(extract_delta(info, fit), error = function(e) { curve_error <<- conditionMessage(e); rep(NA_real_, length(lag_values)) })
    ok <- all(is.finite(prediction)) && all(is.finite(curve))
    list(
      status = if (ok) "success" else "failed",
      prediction = prediction,
      curve = curve,
      cumulative = if (all(is.finite(curve))) mean(curve) else NA_real_,
      singular = isSingular(fit, tol = 1e-4),
      message = paste(c(prediction_error, curve_error)[nzchar(c(prediction_error, curve_error))], collapse = " | ")
    )
  }

  all_idx <- seq_len(nrow(d))
  full_records <- list()
  validation_records <- list()
  prediction_records <- list()

  add_validation <- function(type, heldout, info, structure, result, train_idx, test_idx, include_random_effects) {
    observed <- d$y_log[test_idx]
    predicted <- result$prediction
    valid <- result$status == "success" && length(predicted) == length(observed) && all(is.finite(predicted))
    error <- if (valid) predicted - observed else rep(NA_real_, length(observed))
    validation_records[[length(validation_records) + 1L]] <<- data.table(
      structure = structure$id,
      validation_type = type,
      heldout = as.character(heldout),
      model = info$id,
      fit_status = result$status,
      singular = result$singular,
      n_train = length(train_idx),
      n_test = length(test_idx),
      cumulative_delta = result$cumulative,
      difference_vs_primary = result$cumulative - source_dlm_cumulative,
      rmse_log = if (valid) sqrt(mean(error^2)) else NA_real_,
      mae_log = if (valid) mean(abs(error)) else NA_real_,
      bias_log = if (valid) mean(error) else NA_real_,
      correlation = if (valid) suppressWarnings(cor(predicted, observed)) else NA_real_,
      message = result$message
    )
    held <- d[test_idx]
    prediction_records[[length(prediction_records) + 1L]] <<- data.table(
      structure = structure$id,
      validation_type = type,
      heldout = as.character(heldout),
      model = info$id,
      location_id = held$location_id,
      region_id = held$region_id,
      year = held$year,
      cause_ac = held$cause_ac,
      observed_log_rate = observed,
      predicted_log_rate = if (length(predicted) == length(observed)) predicted else rep(NA_real_, length(observed)),
      prediction_random_effects = include_random_effects
    )
  }

  for (structure in structures) {
    for (info in models) {
      result <- fit_one(info, structure, all_idx, all_idx, TRUE)
      full_records[[length(full_records) + 1L]] <- cbind(
        data.table(structure = structure$id, model = info$id, fit_status = result$status, singular = result$singular, cumulative_delta = result$cumulative, message = result$message),
        as.data.table(as.list(setNames(result$curve, paste0("lag_", lag_values))))
      )
    }
    for (spec in temporal_holdouts) {
      train_idx <- which(d$year <= spec$train_end)
      test_idx <- which(d$year >= spec$test_start & d$year <= spec$test_end)
      for (info in models) add_validation("temporal", spec$id, info, structure, fit_one(info, structure, train_idx, test_idx, TRUE), train_idx, test_idx, TRUE)
    }
  }

  region_ids <- sort(unique(d$region_id))
  if (length(region_ids) != expected_regions || anyNA(region_ids)) stop("Unexpected reference-region count")
  for (structure in structures) {
    for (region_id in region_ids) {
      train_idx <- which(d$region_id != region_id)
      test_idx <- which(d$region_id == region_id)
      for (info in models) add_validation("region", region_id, info, structure, fit_one(info, structure, train_idx, test_idx, FALSE), train_idx, test_idx, FALSE)
    }
  }

  full_fit <- rbindlist(full_records, fill = TRUE)
  validation <- rbindlist(validation_records, fill = TRUE)
  predictions <- rbindlist(prediction_records, fill = TRUE)
  setorder(full_fit, structure, model)
  setorder(validation, structure, validation_type, heldout, model)
  setorder(predictions, structure, validation_type, heldout, model, location_id, year, cause_ac)

  structure_summary <- validation[, .(
    validation_fit_count = .N,
    successful_fit_count = sum(fit_status == "success"),
    singular_fit_count = sum(singular %in% TRUE, na.rm = TRUE),
    all_cumulative_positive = all(cumulative_delta > 0),
    mean_cumulative_delta = mean(cumulative_delta),
    mean_rmse_log = mean(rmse_log)
  ), by = structure]
  temporal_summary <- validation[validation_type == "temporal", .(mean_temporal_rmse_log = mean(rmse_log)), by = structure]
  region_summary <- validation[validation_type == "region", .(mean_region_rmse_log = mean(rmse_log)), by = structure]
  structure_summary <- merge(structure_summary, temporal_summary, by = "structure", all.x = TRUE)
  structure_summary <- merge(structure_summary, region_summary, by = "structure", all.x = TRUE)

  list(full_fit = full_fit, validation = validation, predictions = predictions, structure_summary = structure_summary)
}
