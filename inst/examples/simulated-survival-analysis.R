# Detailed tvcQGComp implementation for the bundled survival dataset
#
# The script estimates three marginal main-effect measures for a simultaneous
# one-quantile increase in all mixture components:
#   HR: pooled person-time marginal hazard ratio
#   RR: end-of-follow-up marginal risk ratio
#   RD: end-of-follow-up marginal risk difference
#
# The point estimates come from the non-bootstrap fit. When the bootstrap is
# enabled, percentile confidence intervals come from bootstrap meta-model
# coefficients transformed to the corresponding effect scale.

# Install the package once from GitHub, if needed:
# install.packages("remotes")
# remotes::install_github("JuwelRana19/tvcQGComp")

library(tvcQGComp)

# ---------------------------------------------------------------------------
# 1. Analysis controls
# ---------------------------------------------------------------------------

# Defaults are suitable for a final analysis. Override them before sourcing:
# Sys.setenv(
#   TVCQGCOMP_MC_SIZE = "1000",
#   TVCQGCOMP_RUN_BOOTSTRAP = "false"
# )
mc_size <- as.integer(Sys.getenv("TVCQGCOMP_MC_SIZE", unset = "10000"))
n_boot <- as.integer(Sys.getenv("TVCQGCOMP_N_BOOT", unset = "200"))
n_workers <- as.integer(Sys.getenv("TVCQGCOMP_N_WORKERS", unset = "4"))
batch_size <- as.integer(
  Sys.getenv("TVCQGCOMP_BATCH_SIZE", unset = as.character(n_workers))
)
run_bootstrap <- identical(
  tolower(Sys.getenv("TVCQGCOMP_RUN_BOOTSTRAP", unset = "false")),
  "true"
)

if (is.na(mc_size) || mc_size < 1L) stop("mc_size must be a positive integer.")
if (is.na(n_boot) || n_boot < 1L) stop("n_boot must be a positive integer.")
if (is.na(n_workers) || n_workers < 1L) {
  stop("n_workers must be a positive integer.")
}
if (is.na(batch_size) || batch_size < 1L) {
  stop("batch_size must be a positive integer.")
}

output_dir <- file.path("tvcQGComp_toy_data_analysis", "detailed_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 2. Load and verify the example data
# ---------------------------------------------------------------------------

data("toy_data", package = "tvcQGComp")
dat <- toy_data

# TimeOut is stored as a factor in the source dataset. The models use numeric
# linear and quadratic time terms.
dat$TimeOut <- as.integer(as.character(dat$TimeOut))

stopifnot(
  nrow(dat) == 5584L,
  length(unique(dat$UniqID)) == 500L,
  sum(dat$status) == 100L,
  all(dat$TimeOut >= dat$TimeInn)
)

# Exposure variables are continuous in the input data. tvcQGComp converts
# each component to quartiles internally before fitting the models.
exposures <- c("BC", "NIT", "SO4", "NH4", "OM")
stopifnot(all(vapply(
  dat[exposures],
  function(values) {
    is.numeric(values) &&
      all(is.finite(values)) &&
      length(unique(values)) > 100L &&
      any(abs(values - round(values)) > 1e-8)
  },
  logical(1)
)))
stopifnot(!any(grepl("_lag[0-9]+$", names(dat))))

# These covariates are simulated forward using the explicit models below.
active_tvc <- c(
  "Urban_form",
  "CanadianRegion",
  "CSize",
  "deprivation",
  "dependency",
  "instability",
  "ethnicconcentration"
)

baseline_covariates <- c(
  "age",
  "income_inadequacy",
  "sex",
  "visible_minority",
  "IndigenousIdentity",
  "marsth",
  "Education",
  "employment",
  "Occupation"
)

factor_variables <- c(
  "income_inadequacy",
  active_tvc,
  baseline_covariates
)

# ---------------------------------------------------------------------------
# 3. Explicit first-stage model formulas
# ---------------------------------------------------------------------------

# Pooled person-time outcome model for the interval event probability.
outcome_formula <- status ~
  BC + NIT + SO4 + NH4 + OM +
  TimeOut + I(TimeOut^2) +
  sex + visible_minority + IndigenousIdentity +
  marsth + Education + employment + Occupation +
  age + income_inadequacy +
  Urban_form + CanadianRegion + CSize +
  deprivation + dependency + instability + ethnicconcentration

# Each active time-varying covariate is predicted from time, baseline
# covariates, and lagged values of the other active TVCs.
tvc_formulas <- list(
  Urban_form = Urban_form ~
    TimeOut + I(TimeOut^2) +
    sex + visible_minority + IndigenousIdentity +
    marsth + Education + employment + Occupation +
    age + income_inadequacy +
    CanadianRegion_lag1 + CSize_lag1 +
    deprivation_lag1 + dependency_lag1 +
    instability_lag1 + ethnicconcentration_lag1,

  CanadianRegion = CanadianRegion ~
    TimeOut + I(TimeOut^2) +
    sex + visible_minority + IndigenousIdentity +
    marsth + Education + employment + Occupation +
    age + income_inadequacy +
    Urban_form_lag1 + CSize_lag1 +
    deprivation_lag1 + dependency_lag1 +
    instability_lag1 + ethnicconcentration_lag1,

  CSize = CSize ~
    TimeOut + I(TimeOut^2) +
    sex + visible_minority + IndigenousIdentity +
    marsth + Education + employment + Occupation +
    age + income_inadequacy +
    Urban_form_lag1 + CanadianRegion_lag1 +
    deprivation_lag1 + dependency_lag1 +
    instability_lag1 + ethnicconcentration_lag1,

  deprivation = deprivation ~
    TimeOut + I(TimeOut^2) +
    sex + visible_minority + IndigenousIdentity +
    marsth + Education + employment + Occupation +
    age + income_inadequacy +
    Urban_form_lag1 + CanadianRegion_lag1 + CSize_lag1 +
    dependency_lag1 + instability_lag1 +
    ethnicconcentration_lag1,

  dependency = dependency ~
    TimeOut + I(TimeOut^2) +
    sex + visible_minority + IndigenousIdentity +
    marsth + Education + employment + Occupation +
    age + income_inadequacy +
    Urban_form_lag1 + CanadianRegion_lag1 + CSize_lag1 +
    deprivation_lag1 + instability_lag1 +
    ethnicconcentration_lag1,

  instability = instability ~
    TimeOut + I(TimeOut^2) +
    sex + visible_minority + IndigenousIdentity +
    marsth + Education + employment + Occupation +
    age + income_inadequacy +
    Urban_form_lag1 + CanadianRegion_lag1 + CSize_lag1 +
    deprivation_lag1 + dependency_lag1 +
    ethnicconcentration_lag1,

  ethnicconcentration = ethnicconcentration ~
    TimeOut + I(TimeOut^2) +
    sex + visible_minority + IndigenousIdentity +
    marsth + Education + employment + Occupation +
    age + income_inadequacy +
    Urban_form_lag1 + CanadianRegion_lag1 + CSize_lag1 +
    deprivation_lag1 + dependency_lag1 + instability_lag1
)

# ---------------------------------------------------------------------------
# 4. Build and validate the tvcQGComp configuration
# ---------------------------------------------------------------------------

meta_targets <- c("HR", "RR", "RD")

cfg <- make_tvcqgcomp_config(
  id = "UniqID",
  time_in = "TimeInn",
  time_out = "TimeOut",
  outcome = "status",
  outcome_type = "survival",
  exposures = exposures,
  time_varying_covariates = active_tvc,
  time_fixed_covariates = baseline_covariates,
  factor_vars = factor_variables,
  outcome_formula = outcome_formula,
  tvc_formulas = tvc_formulas,
  covtypes = c(
    Urban_form = "categorical",
    CanadianRegion = "categorical",
    CSize = "categorical",
    deprivation = "categorical",
    dependency = "categorical",
    instability = "categorical",
    ethnicconcentration = "categorical"
  ),
  q = 4,
  natural_course = "none",
  auto_history = TRUE,
  baselags = TRUE,
  meta_target = meta_targets,
  categorical_sim = "stochastic",
  exposome = "quantized"
)

validate_tvcqgcomp_data(dat, cfg)

# Q0, Q1, Q2, and Q3 assign every exposure to the same quantized level.
intervention_scenarios <- generate_intervention_scenarios(
  config = cfg,
  levels = 0:3,
  prefix = "Q"
)

# ---------------------------------------------------------------------------
# 5. Explicit second-stage marginal meta-models
# ---------------------------------------------------------------------------

# HR uses the pooled person-time hazard outcome and adjusts for follow-up time.
# RR and RD use cumulative risk at the final follow-up time.
meta_formulas <- list(
  HR = hazard ~ joint_effect + factor(TimeInn),
  RR = cumrisk ~ joint_effect,
  RD = cumrisk ~ joint_effect
)

meta_families <- list(
  HR = stats::quasibinomial(link = "logit"),
  RR = stats::quasibinomial(link = "log"),
  RD = stats::quasibinomial(link = "identity")
)

# ---------------------------------------------------------------------------
# 6. Fit non-bootstrap point estimates
# ---------------------------------------------------------------------------

fit_started <- Sys.time()

fit <- tvcQGComp_survival(
  data = dat,
  config = cfg,
  mc_size = mc_size,
  intervention_scenarios = intervention_scenarios,
  natural_course = "none",
  seed = 1234L,
  replace_mc = TRUE,
  meta_formula = meta_formulas,
  meta_family = meta_families,
  meta_target = meta_targets,
  fit_meta = TRUE,
  fit_exposure_models = FALSE,
  exposure_scale = "quantized",
  verbose = TRUE
)

fit_finished <- Sys.time()

# Model-based effects per simultaneous one-quantile increase in the mixture:
# HR and RR are exp(beta_joint_effect); RD is beta_joint_effect.
point_effect_summary <- as.data.frame(fit$meta_effect_summary)

# Retain every fitted meta-model coefficient for auditing.
meta_model_coefficients <- do.call(
  rbind,
  lapply(names(fit$meta_models), function(target) {
    coefficient <- stats::coef(fit$meta_models[[target]])
    data.frame(
      meta_target = target,
      term = names(coefficient),
      link_scale_estimate = unname(coefficient),
      row.names = NULL
    )
  })
)

cat("\nModel-based marginal effect estimates\n")
print(point_effect_summary)

cat(
  "\nNon-bootstrap elapsed minutes:",
  round(as.numeric(difftime(fit_finished, fit_started, units = "mins")), 2),
  "\n"
)

# ---------------------------------------------------------------------------
# 7. Diagnostics, trajectories, and non-bootstrap outputs
# ---------------------------------------------------------------------------

diagnostics <- diagnose_tvcqgcomp_run(fit)
risk_trajectory <- as.data.frame(fit$risk_trajectory)

utils::write.csv(
  point_effect_summary,
  file.path(output_dir, "meta_effect_summary_point_estimates.csv"),
  row.names = FALSE
)
utils::write.csv(
  meta_model_coefficients,
  file.path(output_dir, "meta_model_coefficients.csv"),
  row.names = FALSE
)
utils::write.csv(
  risk_trajectory,
  file.path(output_dir, "risk_trajectory.csv"),
  row.names = FALSE
)

# Uncompressed storage avoids a long xz-compression delay after model fitting.
saveRDS(fit, file.path(output_dir, "tvcQGComp_fit.rds"), compress = FALSE)
saveRDS(
  diagnostics,
  file.path(output_dir, "tvcQGComp_diagnostics.rds"),
  compress = FALSE
)

# The primary figure is cumulative incidence by quantile.
grDevices::png(
  file.path(output_dir, "cumulative_risk_trajectories.png"),
  width = 2400,
  height = 1500,
  res = 300
)
plot_cumulative_risk_trajectory(fit, include_natural = FALSE, include_average = TRUE)
grDevices::dev.off()

# Survival is an optional separate figure.
grDevices::png(
  file.path(output_dir, "survival_trajectories.png"),
  width = 2400,
  height = 1500,
  res = 300
)
plot_survival_trajectory(fit, include_natural = FALSE, include_average = TRUE)
grDevices::dev.off()

# ---------------------------------------------------------------------------
# 8. Cluster bootstrap confidence intervals
# ---------------------------------------------------------------------------

if (run_bootstrap) {
  checkpoint_file <- file.path(output_dir, "bootstrap_checkpoint.rds")
  bootstrap_started <- Sys.time()

  boot <- tvcQGComp_survival_boot(
    data = dat,
    config = cfg,
    n_boot = n_boot,
    mc_size = mc_size,
    intervention_scenarios = intervention_scenarios,
    natural_course = "none",
    seed = 1234L,
    replace_mc = TRUE,
    meta_formula = meta_formulas,
    meta_family = meta_families,
    meta_target = meta_targets,
    fit_meta = TRUE,
    fit_exposure_models = FALSE,
    exposure_scale = "quantized",
    parallel = n_workers > 1L,
    n_workers = n_workers,
    batch_size = batch_size,
    checkpoint_file = checkpoint_file,
    resume_from_checkpoint = TRUE,
    verbose = TRUE,
    stop_on_error = FALSE,
    keep_boot_objects = FALSE,
    keep_boot_curves = TRUE
  )

  bootstrap_finished <- Sys.time()

  boot_results <- if (
    is.list(boot) &&
      !inherits(boot, "data.frame") &&
      !is.null(boot$results)
  ) {
    as.data.frame(boot$results)
  } else {
    as.data.frame(boot)
  }

  successful_boot <- boot_results[boot_results$status == "success", , drop = FALSE]
  if (!nrow(successful_boot)) {
    stop("No bootstrap replicate completed successfully.")
  }

  # Point estimates remain from the original non-bootstrap fit. Confidence
  # limits are percentile intervals from bootstrap meta-model coefficients.
  bootstrap_effect_summary <- do.call(
    rbind,
    lapply(meta_targets, function(target) {
      coefficient_column <- paste0("coef_", target, "_joint_effect")
      if (!(coefficient_column %in% names(successful_boot))) {
        stop("Bootstrap output is missing ", coefficient_column, ".")
      }

      bootstrap_link_estimates <- successful_boot[[coefficient_column]]
      bootstrap_effect_estimates <- if (target %in% c("HR", "RR")) {
        exp(bootstrap_link_estimates)
      } else {
        bootstrap_link_estimates
      }

      point_row <- point_effect_summary[
        point_effect_summary$meta_target == target,
        ,
        drop = FALSE
      ]
      percentile_ci <- stats::quantile(
        bootstrap_effect_estimates,
        probs = c(0.025, 0.975),
        na.rm = TRUE,
        names = FALSE
      )

      data.frame(
        meta_target = target,
        point_estimate = point_row$estimate,
        bootstrap_se = stats::sd(
          bootstrap_effect_estimates,
          na.rm = TRUE
        ),
        ci_lower_95 = percentile_ci[[1L]],
        ci_upper_95 = percentile_ci[[2L]],
        n_success = sum(is.finite(bootstrap_effect_estimates)),
        n_requested = n_boot,
        row.names = NULL
      )
    })
  )

  cat("\nPoint estimates with bootstrap percentile confidence intervals\n")
  print(bootstrap_effect_summary)
  cat(
    "\nBootstrap elapsed minutes:",
    round(
      as.numeric(
        difftime(bootstrap_finished, bootstrap_started, units = "mins")
      ),
      2
    ),
    "\n"
  )

  utils::write.csv(
    boot_results,
    file.path(output_dir, "bootstrap_meta_model_coefficients.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    bootstrap_effect_summary,
    file.path(output_dir, "meta_effect_summary_bootstrap_ci.csv"),
    row.names = FALSE
  )
  saveRDS(boot, file.path(output_dir, "tvcQGComp_bootstrap.rds"))
}

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "session_info.txt")
)

cat("\nOutputs saved to:", normalizePath(output_dir), "\n")
