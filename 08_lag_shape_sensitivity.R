suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(lme4)
})

# Registered lag-shape and deterministic lower/upper sensitivity scenarios.
# Requires 00_gbd_analysis_core.R. Primary curve and cumulative estimate are
# passed directly from the primary analysis result.

run_lag_shape_sensitivity <- function(outcome, exposure, sdi, units, primary_curve, primary_cumulative) {
  prepared <- prepare_gbd_dlm_data(outcome, exposure, sdi, units)
  y <- prepared$outcome
  analysis_years <- prepared$analysis_years
  lag_values <- prepared$lag_values
  lag_cols <- prepared$lag_cols
  if (length(primary_curve) != length(lag_values) || !all(is.finite(primary_curve)) || !is.finite(primary_cumulative)) {
    stop("Primary DLM reference must contain 11 finite lag estimates and one finite cumulative estimate")
  }

  scenarios <- list(
    list(id = "baseline_lme4", outcome_field = "val", lag_df = 4L, lag_knots = c(3, 7), calendar = "natural_spline_df3", sdi = "residualized_calendar"),
    list(id = "simpler_lag_df3", outcome_field = "val", lag_df = 3L, lag_knots = 5, calendar = "natural_spline_df3", sdi = "residualized_calendar"),
    list(id = "no_sdi", outcome_field = "val", lag_df = 4L, lag_knots = c(3, 7), calendar = "natural_spline_df3", sdi = "none"),
    list(id = "linear_calendar", outcome_field = "val", lag_df = 4L, lag_knots = c(3, 7), calendar = "linear", sdi = "residualized_calendar"),
    list(id = "ui_lower", outcome_field = "lower", lag_df = 4L, lag_knots = c(3, 7), calendar = "natural_spline_df3", sdi = "residualized_calendar"),
    list(id = "ui_upper", outcome_field = "upper", lag_df = 4L, lag_knots = c(3, 7), calendar = "natural_spline_df3", sdi = "residualized_calendar")
  )

  fit_scenario <- function(sc) {
    response <- as.numeric(y[[sc$outcome_field]])
    if (any(!is.finite(response)) || any(response <= 0)) stop("Scenario outcome values must be positive and finite")
    lag_spline <- ns(lag_values, df = sc$lag_df, knots = sc$lag_knots, Boundary.knots = c(0, 10), intercept = TRUE)
    lag_raw <- as.matrix(y[, ..lag_cols]) %*% lag_spline
    lag_center <- colMeans(lag_raw)
    lag_qr <- qr(sweep(lag_raw, 2L, lag_center, "-"), LAPACK = FALSE)
    lag_scale <- sqrt(nrow(lag_raw) - 1)
    lag_orth <- qr.Q(lag_qr) * lag_scale
    lag_transform <- qr.R(lag_qr) / lag_scale
    lag_names <- paste0("lag_q", seq_len(ncol(lag_orth)))
    interaction_names <- paste0("cause_lag_q", seq_len(ncol(lag_orth)))

    if (sc$calendar == "natural_spline_df3") {
      time_features <- ns(analysis_years, df = 3, Boundary.knots = c(min(analysis_years), max(analysis_years)))[y$time_idx, , drop = FALSE]
    } else {
      time_features <- matrix(y$time_std, ncol = 1)
    }
    time_names <- paste0("time_q", seq_len(ncol(time_features)))
    sdi_conditional <- if (sc$sdi == "none") NULL else as.numeric(residuals(lm(y$sdi_std ~ time_features)))

    d <- data.table(
      y_log = log(response),
      cause_ac = y$cause_ac,
      time_std = y$time_std,
      location_factor = y$location_factor
    )
    for (k in seq_len(ncol(lag_orth))) d[[lag_names[k]]] <- lag_orth[, k]
    for (k in seq_len(ncol(lag_orth))) d[[interaction_names[k]]] <- d$cause_ac * d[[lag_names[k]]]
    for (j in seq_len(ncol(time_features))) d[[time_names[j]]] <- time_features[, j]
    fixed_terms <- c("cause_ac", lag_names, interaction_names)
    if (!is.null(sdi_conditional)) {
      d[, sdi_conditional := sdi_conditional]
      fixed_terms <- c(fixed_terms, "sdi_conditional")
    }
    fixed_terms <- c(fixed_terms, time_names)
    formula <- as.formula(paste("y_log ~", paste(fixed_terms, collapse = " + "), "+ (1 + time_std | location_factor)"))
    fitted <- fit_lmer_with_fallback(formula, d)
    if (is.null(fitted$fit)) return(list(status = "failed", curve = rep(NA_real_, length(lag_values)), cumulative = NA_real_, singular = NA, message = fitted$message))
    fit <- fitted$fit
    fixed <- fixef(fit)
    if (!all(interaction_names %in% names(fixed))) return(list(status = "failed", curve = rep(NA_real_, length(lag_values)), cumulative = NA_real_, singular = isSingular(fit, tol = 1e-4), message = "Missing interaction coefficients"))
    gamma <- as.numeric(fixed[interaction_names])
    delta_original <- as.numeric(solve(lag_transform, gamma))
    curve <- as.numeric(lag_spline %*% delta_original)
    list(status = "success", curve = curve, cumulative = mean(curve), singular = isSingular(fit, tol = 1e-4), message = "")
  }

  summary_rows <- list()
  delta_rows <- list()
  for (sc in scenarios) {
    result <- fit_scenario(sc)
    curve <- result$curve
    ok <- result$status == "success" && all(is.finite(curve))
    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      scenario = sc$id,
      fit_status = result$status,
      outcome_field = sc$outcome_field,
      lag_df = sc$lag_df,
      calendar = sc$calendar,
      sdi = sc$sdi,
      cumulative_delta = result$cumulative,
      cumulative_delta_difference_vs_primary = if (ok) result$cumulative - primary_cumulative else NA_real_,
      delta_curve_rmse_vs_primary = if (ok) sqrt(mean((curve - primary_curve)^2)) else NA_real_,
      delta_curve_sign_consistency = if (ok) mean(sign(curve) == sign(primary_curve)) else NA_real_,
      singular = result$singular,
      message = result$message
    )
    delta_rows[[length(delta_rows) + 1L]] <- data.table(
      scenario = sc$id,
      lag = lag_values,
      delta_sensitivity = curve,
      delta_primary = primary_curve,
      difference = curve - primary_curve,
      sign_consistent = sign(curve) == sign(primary_curve)
    )
  }
  summary <- rbindlist(summary_rows, fill = TRUE)
  delta <- rbindlist(delta_rows, fill = TRUE)
  if (!all(summary$fit_status == "success")) stop("One or more registered lag-shape sensitivity scenarios failed")
  list(summary = summary, delta = delta)
}
