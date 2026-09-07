suppressPackageStartupMessages({
  library(data.table)
  library(rjags)
  library(coda)
  library(lme4)
})

# Primary 0-10-year cumulative AC-versus-OC association.
# Requires prepare_gbd_dlm_data() from 00_gbd_analysis_core.R to be loaded.
# All inputs are in-memory data frames; no project paths or output files are used.

run_primary_cumulative_model <- function(
  outcome,
  exposure,
  sdi,
  units,
  bootstrap_replicates = 300L,
  seed = 20260903L,
  minimum_successful_bootstrap = 285L
) {
  prepared <- prepare_gbd_dlm_data(outcome, exposure, sdi, units)
  y <- prepared$outcome
  lag_values <- prepared$lag_values
  lag_features <- prepared$lag_features
  time_features <- prepared$time_features
  lag_transform <- prepared$lag_transform
  lag_spline <- prepared$lag_spline

  location_levels <- sort(unique(y$location_id))
  location_map <- setNames(seq_along(location_levels), location_levels)
  y[, loc_idx := unname(location_map[as.character(location_id)])]

  model_text <- paste(c(
    "model {",
    "  for (n in 1:N) {",
    "    y[n] ~ dnorm(mu[n], tau_e)",
    "    mu[n] <- alpha + alpha_ac * cause_ac[n] + inprod(beta_orth[1:K], lag_features[n,1:K]) + cause_ac[n] * inprod(delta_orth[1:K], lag_features[n,1:K]) + sdi_beta * sdi_std[n] + inprod(time_beta[1:J], time_features[n,1:J]) + b[loc_idx[n]] + q[loc_idx[n]] * time_std[n]",
    "  }",
    "  alpha ~ dnorm(0, 0.01)",
    "  alpha_ac ~ dnorm(0, 1)",
    "  sdi_beta ~ dnorm(0, 1)",
    "  for (k in 1:K) { beta_orth[k] ~ dnorm(0, 1); delta_orth[k] ~ dnorm(0, 1) }",
    "  for (j in 1:J) { time_beta[j] ~ dnorm(0, 1) }",
    "  sigma_e ~ dunif(0, 2); tau_e <- pow(sigma_e, -2)",
    "  sigma_b ~ dunif(0, 2); tau_b <- pow(sigma_b, -2)",
    "  sigma_q ~ dunif(0, 1); tau_q <- pow(sigma_q, -2)",
    "  for (i in 1:n_loc) { b[i] ~ dnorm(0, tau_b); q[i] ~ dnorm(0, tau_q) }",
    "}"
  ), collapse = "\n")
  model_connection <- textConnection(model_text)
  on.exit(close(model_connection), add = TRUE)

  jags_data <- list(
    N = nrow(y),
    K = ncol(lag_features),
    J = ncol(time_features),
    n_loc = length(location_levels),
    y = y$y_log,
    cause_ac = y$cause_ac,
    loc_idx = y$loc_idx,
    time_std = y$time_std,
    sdi_std = prepared$sdi_conditional,
    lag_features = lag_features,
    time_features = time_features
  )

  set.seed(seed)
  jm <- jags.model(model_connection, data = jags_data, n.chains = 3L, quiet = TRUE)
  update(jm, n.iter = 1000L, progress.bar = "none")
  monitored <- c(
    "alpha", "alpha_ac", "sdi_beta",
    paste0("beta_orth[", seq_len(ncol(lag_features)), "]"),
    paste0("delta_orth[", seq_len(ncol(lag_features)), "]"),
    paste0("time_beta[", seq_len(ncol(time_features)), "]"),
    "sigma_e", "sigma_b", "sigma_q"
  )
  samples <- coda.samples(jm, variable.names = monitored, n.iter = 1800L, thin = 2L, progress.bar = "none")

  beta_orth_names <- paste0("beta_orth[", seq_len(ncol(lag_features)), "]")
  delta_orth_names <- paste0("delta_orth[", seq_len(ncol(lag_features)), "]")
  scalar_names <- c("alpha", "alpha_ac", "sdi_beta", paste0("time_beta[", seq_len(ncol(time_features)), "]"), "sigma_e", "sigma_b", "sigma_q")
  posterior <- mcmc.list(lapply(samples, function(chain) {
    m <- as.matrix(chain)
    beta_original <- m[, beta_orth_names, drop = FALSE] %*% t(solve(lag_transform))
    delta_original <- m[, delta_orth_names, drop = FALSE] %*% t(solve(lag_transform))
    colnames(beta_original) <- paste0("beta[", seq_len(ncol(beta_original)), "]")
    colnames(delta_original) <- paste0("delta[", seq_len(ncol(delta_original)), "]")
    as.mcmc(cbind(m[, scalar_names, drop = FALSE], beta_original, delta_original))
  }))

  posterior_matrix <- as.matrix(posterior)
  delta_names <- paste0("delta[", seq_len(ncol(lag_features)), "]")
  delta_draws <- posterior_matrix[, delta_names, drop = FALSE]
  delta_by_lag <- delta_draws %*% t(lag_spline)
  delta_mcmc <- mcmc.list(lapply(posterior, function(chain) {
    z <- as.matrix(chain)[, delta_names, drop = FALSE] %*% t(lag_spline)
    colnames(z) <- paste0("lag_", lag_values)
    as.mcmc(z)
  }))
  cumulative_draws <- rowMeans(delta_by_lag)
  cumulative_mcmc <- mcmc.list(lapply(delta_mcmc, function(chain) {
    as.mcmc(matrix(rowMeans(as.matrix(chain)), ncol = 1, dimnames = list(NULL, "cumulative_0_10")))
  }))

  delta_summary <- rbindlist(lapply(seq_along(lag_values), function(j) {
    z <- delta_by_lag[, j]
    q <- quantile(z, c(0.025, 0.5, 0.975), names = FALSE)
    data.table(lag = lag_values[j], mean = mean(z), lower = q[1], median = q[2], upper = q[3])
  }))
  cumulative_q <- quantile(cumulative_draws, c(0.025, 0.5, 0.975), names = FALSE)
  delta_summary <- rbind(
    delta_summary,
    data.table(lag = NA_integer_, mean = mean(cumulative_draws), lower = cumulative_q[1], median = cumulative_q[2], upper = cumulative_q[3])
  )

  # Country/territory block bootstrap, preserving complete AC/OC trajectories.
  d <- copy(prepared$d)
  lag_names <- prepared$lag_names
  interaction_names <- prepared$interaction_names
  time_names <- prepared$time_names
  fixed_terms <- c("cause_ac", lag_names, interaction_names, "sdi_conditional", time_names)
  bootstrap_formula <- as.formula(paste("y_log ~", paste(fixed_terms, collapse = " + "), "+ (1 + time_std | boot_loc)"))
  unit_rows <- split(d, d$location_id)
  unit_ids <- names(unit_rows)
  set.seed(seed)
  bootstrap_records <- vector("list", bootstrap_replicates)

  for (b in seq_len(bootstrap_replicates)) {
    sampled <- sample(unit_ids, length(unit_ids), replace = TRUE)
    boot <- rbindlist(lapply(seq_along(sampled), function(j) {
      z <- copy(unit_rows[[sampled[j]]])
      z[, boot_loc := paste0(sampled[j], "_", j)]
      z
    }))
    boot[, boot_loc := factor(boot_loc)]
    result <- tryCatch({
      fit <- suppressWarnings(lmer(
        bootstrap_formula,
        data = boot,
        REML = FALSE,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 200000),
          check.rankX = "silent.drop.cols",
          check.conv.singular = "ignore"
        )
      ))
      fixed <- fixef(fit)
      gamma <- as.numeric(fixed[interaction_names])
      if (anyNA(gamma)) stop("Missing bootstrap interaction coefficient")
      delta_original <- as.numeric(solve(lag_transform, gamma))
      curve <- as.numeric(lag_spline %*% delta_original)
      list(status = "success", curve = curve, cumulative = mean(curve), singular = isSingular(fit, tol = 1e-4), message = "")
    }, error = function(e) {
      list(status = "failed", curve = rep(NA_real_, length(lag_values)), cumulative = NA_real_, singular = NA, message = conditionMessage(e))
    })
    bootstrap_records[[b]] <- cbind(
      data.table(replicate = b, status = result$status, cumulative = result$cumulative, singular = result$singular, message = result$message),
      as.data.table(as.list(setNames(result$curve, paste0("lag_", lag_values))))
    )
  }

  bootstrap <- rbindlist(bootstrap_records, fill = TRUE)
  valid_bootstrap <- bootstrap[status == "success" & is.finite(cumulative)]
  if (nrow(valid_bootstrap) < minimum_successful_bootstrap) stop("Bootstrap successful replicate count below the registered minimum")
  boot_curves <- as.matrix(valid_bootstrap[, paste0("lag_", lag_values), with = FALSE])
  bootstrap_summary <- rbindlist(lapply(seq_along(lag_values), function(j) {
    q <- quantile(boot_curves[, j], c(0.025, 0.5, 0.975), names = FALSE)
    data.table(lag = lag_values[j], bootstrap_lower = q[1], bootstrap_median = q[2], bootstrap_upper = q[3])
  }))
  boot_cum_q <- quantile(valid_bootstrap$cumulative, c(0.025, 0.5, 0.975), names = FALSE)
  bootstrap_summary <- rbind(
    bootstrap_summary,
    data.table(lag = NA_integer_, bootstrap_lower = boot_cum_q[1], bootstrap_median = boot_cum_q[2], bootstrap_upper = boot_cum_q[3])
  )

  delta <- merge(delta_summary, bootstrap_summary, by = "lag", all = TRUE, sort = TRUE)
  summary <- data.table(
    parameter = c("posterior_cumulative", "bootstrap_cumulative", "sigma_e", "sigma_b", "sigma_q"),
    mean = c(mean(cumulative_draws), mean(valid_bootstrap$cumulative), mean(posterior_matrix[, "sigma_e"]), mean(posterior_matrix[, "sigma_b"]), mean(posterior_matrix[, "sigma_q"])),
    lower = c(cumulative_q[1], boot_cum_q[1], quantile(posterior_matrix[, "sigma_e"], 0.025), quantile(posterior_matrix[, "sigma_b"], 0.025), quantile(posterior_matrix[, "sigma_q"], 0.025)),
    upper = c(cumulative_q[3], boot_cum_q[3], quantile(posterior_matrix[, "sigma_e"], 0.975), quantile(posterior_matrix[, "sigma_b"], 0.975), quantile(posterior_matrix[, "sigma_q"], 0.975))
  )

  delta_psrf <- gelman.diag(delta_mcmc, autoburnin = FALSE, multivariate = FALSE)$psrf[, "Point est."]
  delta_ess <- effectiveSize(delta_mcmc)
  cumulative_psrf <- gelman.diag(cumulative_mcmc, autoburnin = FALSE, multivariate = FALSE)$psrf[, "Point est."]
  cumulative_ess <- effectiveSize(cumulative_mcmc)
  diagnostics <- list(
    locations = length(location_levels),
    observations = nrow(y),
    posterior_draws = nrow(posterior_matrix),
    chains = 3L,
    primary_max_rhat = max(c(delta_psrf, cumulative_psrf), na.rm = TRUE),
    primary_min_effective_size = min(c(delta_ess, cumulative_ess), na.rm = TRUE),
    bootstrap_requested = bootstrap_replicates,
    bootstrap_success = nrow(valid_bootstrap),
    bootstrap_singular = sum(valid_bootstrap$singular %in% TRUE, na.rm = TRUE),
    uncertainty_scope = "Posterior and country-block bootstrap conditional on GBD point estimates"
  )

  list(delta = delta, bootstrap = bootstrap, summary = summary, diagnostics = diagnostics, posterior = posterior)
}
