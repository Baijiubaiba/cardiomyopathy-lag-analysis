suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(lme4)
})

# Same-year WHO alcohol-consumption exposure-source triangulation.
# Outcomes and SDI remain GBD-derived; this is not independent outcome validation.
# All input tables are supplied in memory.

run_who_exposure_triangulation <- function(
  outcome,
  sdi,
  who,
  units,
  country_crosswalk,
  bootstrap_replicates = 300L,
  seed = 20260907L,
  minimum_successful_bootstrap = 285L
) {
  units <- assert_primary_units(units)
  primary_ids <- units$location_id
  crosswalk <- as.data.table(copy(country_crosswalk))
  who <- as.data.table(copy(who))
  outcome <- as.data.table(copy(outcome))
  sdi <- as.data.table(copy(sdi))

  required_crosswalk <- c("gbd_location_id", "gbd_location_name", "iso3", "who_code_present")
  required_who <- c("IND_CODE", "DIM_TIME", "DIM_TIME_TYPE", "DIM_GEO_CODE_M49", "DIM_GEO_CODE_TYPE", "DIM_PUBLISH_STATE_CODE", "GEO_NAME_SHORT", "DIM_SEX", "RATE_PER_CAPITA_N")
  if (!all(required_crosswalk %in% names(crosswalk))) stop("country_crosswalk is missing required fields")
  if (!all(required_who %in% names(who))) stop("WHO exposure table is missing required fields")

  crosswalk[, gbd_location_id := as.integer(gbd_location_id)]
  crosswalk <- crosswalk[gbd_location_id %in% primary_ids]
  if (nrow(crosswalk) != 204L || uniqueN(crosswalk$gbd_location_id) != 204L) stop("Crosswalk must cover the 204 registered units")

  who[, c("IND_CODE", "DIM_TIME_TYPE", "DIM_GEO_CODE_TYPE", "DIM_PUBLISH_STATE_CODE", "GEO_NAME_SHORT", "DIM_SEX") :=
    lapply(.SD, as.character), .SDcols = c("IND_CODE", "DIM_TIME_TYPE", "DIM_GEO_CODE_TYPE", "DIM_PUBLISH_STATE_CODE", "GEO_NAME_SHORT", "DIM_SEX")]
  who[, c("DIM_TIME", "DIM_GEO_CODE_M49") := list(as.integer(DIM_TIME), as.integer(DIM_GEO_CODE_M49))]
  who[, RATE_PER_CAPITA_N := as.numeric(RATE_PER_CAPITA_N)]
  who <- who[
    IND_CODE == "SA_0000001688" & DIM_TIME_TYPE == "YEAR" & DIM_GEO_CODE_TYPE == "COUNTRY" &
      DIM_PUBLISH_STATE_CODE == "PUBLISHED" & DIM_SEX == "TOTAL" & DIM_TIME %in% 2000:2022
  ]
  if (anyDuplicated(who[, .(DIM_GEO_CODE_M49, DIM_TIME)])) stop("WHO exposure contains duplicate country-year keys")
  if (any(!is.finite(who$RATE_PER_CAPITA_N))) stop("WHO exposure contains non-finite values")

  normalize_name <- function(x) {
    y <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
    tolower(gsub("[^a-z0-9]", "", tolower(y)))
  }
  who[, name_key := normalize_name(GEO_NAME_SHORT)]
  who[name_key == "netherlandskingdomofthe", name_key := "netherlands"]
  who[name_key == "unitedkingdomofgreatbritainandnorthernireland", name_key := "unitedkingdom"]
  who_counts <- who[, .(who_n_years = uniqueN(DIM_TIME), who_name = first(GEO_NAME_SHORT), who_m49 = first(DIM_GEO_CODE_M49)), by = name_key]
  crosswalk[, name_key := normalize_name(gbd_location_name)]
  crosswalk <- merge(crosswalk, who_counts, by = "name_key", all.x = TRUE, sort = FALSE)
  crosswalk[, included := !is.na(who_m49) & who_n_years == 23L]
  if (crosswalk[, sum(included)] != 188L) stop("Expected 188 units with complete WHO exposure coverage")

  external_ids <- crosswalk[included == TRUE, gbd_location_id]
  coverage <- crosswalk[, .(location_id = gbd_location_id, gbd_location_name, iso3, who_name, who_m49, who_n_years, included)]

  outcome[, c("location_id", "year", "cause_id", "measure_id", "metric_id", "age_id", "sex_id") := lapply(.SD, as.integer), .SDcols = c("location_id", "year", "cause_id", "measure_id", "metric_id", "age_id", "sex_id")]
  sdi[, c("location_id", "year_id", "age_group_id", "sex_id") := lapply(.SD, as.integer), .SDcols = c("location_id", "year_id", "age_group_id", "sex_id")]
  analysis_years <- 2000:2022
  outcome <- outcome[cause_id %in% c(938L, 944L) & measure_id == 2L & metric_id == 3L & age_id == 22L & sex_id == 3L & year %in% analysis_years & location_id %in% external_ids]
  sdi <- sdi[covariate_name_short == "sdi" & age_group_id == 22L & sex_id == 3L & year_id %in% analysis_years & location_id %in% external_ids]
  if (nrow(outcome) != 188L * 23L * 2L || nrow(sdi) != 188L * 23L) stop("Unexpected GBD overlap coverage")
  if (any(!is.finite(outcome$val)) || any(outcome$val <= 0) || any(!is.finite(sdi$mean_value))) stop("Invalid GBD point estimates")

  who_lookup <- merge(
    crosswalk[included == TRUE, .(location_id = gbd_location_id, name_key)],
    who[, .(name_key, year = DIM_TIME, who_apc = RATE_PER_CAPITA_N)],
    by = "name_key", allow.cartesian = FALSE, sort = FALSE
  )
  if (nrow(who_lookup) != 188L * 23L || who_lookup[, .N, by = .(location_id, year)][, any(N != 1L)]) stop("WHO analysis lookup is incomplete or duplicated")
  who_lookup[, name_key := NULL]
  sdi_lookup <- sdi[, .(location_id, year = year_id, sdi = as.numeric(mean_value))]
  outcome <- merge(outcome, who_lookup, by = c("location_id", "year"), all.x = TRUE, sort = FALSE)
  outcome <- merge(outcome, sdi_lookup, by = c("location_id", "year"), all.x = TRUE, sort = FALSE)
  if (anyNA(outcome$who_apc) || anyNA(outcome$sdi)) stop("Incomplete WHO/GBD merged analysis table")

  who_stats <- who_lookup[, .(mean = mean(who_apc), sd = sd(who_apc))]
  sdi_stats <- sdi_lookup[, .(mean = mean(sdi), sd = sd(sdi))]
  if (who_stats$sd <= 0 || sdi_stats$sd <= 0) stop("Standardization scales must be positive")
  outcome[, who_apc_std := (who_apc - who_stats$mean) / who_stats$sd]
  outcome[, sdi_std := (sdi - sdi_stats$mean) / sdi_stats$sd]
  outcome[, time_idx := match(year, analysis_years)]
  outcome[, time_std := 2 * (time_idx - 1) / (length(analysis_years) - 1) - 1]
  outcome[, cause_ac := as.integer(cause_id == 938L)]
  outcome[, y_log := log(as.numeric(val))]
  time_features <- ns(analysis_years, df = 3L, Boundary.knots = c(2000, 2022))[outcome$time_idx, , drop = FALSE]
  sdi_conditional <- as.numeric(residuals(lm(outcome$sdi_std ~ time_features)))

  d <- data.table(
    location_id = outcome$location_id,
    year = outcome$year,
    y_log = outcome$y_log,
    cause_ac = outcome$cause_ac,
    who_apc_std = outcome$who_apc_std,
    sdi_conditional = sdi_conditional,
    time_std = outcome$time_std
  )
  d[, location_factor := factor(as.character(location_id), levels = as.character(sort(unique(location_id))))]
  d[, cause_who_apc := cause_ac * who_apc_std]
  time_names <- paste0("time_q", seq_len(ncol(time_features)))
  for (j in seq_len(ncol(time_features))) d[[time_names[j]]] <- time_features[, j]

  structures <- c("random_intercept_slope", "random_intercept_only")
  fit_one <- function(data, group_name, structure) {
    random_term <- if (structure == "random_intercept_slope") paste0("(1 + time_std |", group_name, ")") else paste0("(1 |", group_name, ")")
    formula <- as.formula(paste("y_log ~ cause_ac + who_apc_std + cause_who_apc + sdi_conditional +", paste(time_names, collapse = " + "), "+", random_term))
    fitted <- fit_lmer_with_fallback(formula, data)
    if (is.null(fitted$fit)) return(list(status = "failed", estimate = NA_real_, ratio = NA_real_, se = NA_real_, singular = NA, message = fitted$message))
    fit <- fitted$fit
    fixed <- fixef(fit)
    if (!("cause_who_apc" %in% names(fixed)) || !is.finite(fixed[["cause_who_apc"]])) return(list(status = "failed", estimate = NA_real_, ratio = NA_real_, se = NA_real_, singular = isSingular(fit, tol = 1e-4), message = "Missing cause-by-WHO interaction"))
    estimate <- as.numeric(fixed[["cause_who_apc"]])
    vc <- tryCatch(vcov(fit), error = function(e) NULL)
    se <- if (!is.null(vc) && "cause_who_apc" %in% rownames(vc)) sqrt(as.numeric(vc["cause_who_apc", "cause_who_apc"])) else NA_real_
    list(status = "success", estimate = estimate, ratio = exp(estimate), se = se, singular = isSingular(fit, tol = 1e-4), message = "")
  }

  full_results <- setNames(lapply(structures, function(s) fit_one(d, "location_factor", s)), structures)
  full <- rbindlist(lapply(structures, function(s) {
    r <- full_results[[s]]
    data.table(structure = s, fit_status = r$status, singular = r$singular, estimate_log_rate_contrast = r$estimate, model_implied_ratio = r$ratio, standard_error = r$se, message = r$message)
  }))

  unit_rows <- split(d, d$location_id)
  unit_ids <- names(unit_rows)
  set.seed(seed)
  records <- vector("list", bootstrap_replicates * length(structures))
  k <- 0L
  for (b in seq_len(bootstrap_replicates)) {
    sampled <- sample(unit_ids, length(unit_ids), replace = TRUE)
    boot <- rbindlist(lapply(seq_along(sampled), function(j) {
      z <- copy(unit_rows[[sampled[j]]])
      z[, boot_loc := paste0(sampled[j], "_", j)]
      z
    }), fill = TRUE)
    boot[, boot_loc := factor(boot_loc)]
    for (s in structures) {
      r <- fit_one(boot, "boot_loc", s)
      k <- k + 1L
      records[[k]] <- data.table(replicate = b, structure = s, status = r$status, estimate_log_rate_contrast = r$estimate, ratio = r$ratio, singular = r$singular, message = r$message)
    }
  }
  bootstrap <- rbindlist(records, fill = TRUE)
  summary <- rbindlist(lapply(structures, function(s) {
    full_result <- full_results[[s]]
    valid <- bootstrap[structure == s & status == "success" & is.finite(estimate_log_rate_contrast)]
    if (nrow(valid) < minimum_successful_bootstrap) stop("WHO bootstrap success below registered minimum")
    q <- quantile(valid$estimate_log_rate_contrast, c(0.025, 0.5, 0.975), names = FALSE)
    ratio_q <- quantile(valid$ratio, c(0.025, 0.5, 0.975), names = FALSE)
    data.table(
      structure = s,
      estimate_log_rate_contrast = full_result$estimate,
      model_implied_ratio = full_result$ratio,
      bootstrap_success = nrow(valid),
      bootstrap_lower = q[1],
      bootstrap_median = q[2],
      bootstrap_upper = q[3],
      bootstrap_ratio_lower = ratio_q[1],
      bootstrap_ratio_median = ratio_q[2],
      bootstrap_ratio_upper = ratio_q[3],
      full_fit_singular = full_result$singular,
      bootstrap_singular = sum(valid$singular %in% TRUE, na.rm = TRUE),
      direction = ifelse(full_result$estimate > 0, "positive", "non_positive")
    )
  }))

  list(coverage = coverage, full = full, bootstrap = bootstrap, summary = summary)
}
