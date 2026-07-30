with_formula_helpers <- function(formula) {
  if (is.null(formula)) {
    return(NULL)
  }
  if (!inherits(formula, "formula")) {
    stop("Expected a formula object.")
  }

  base_env <- environment(formula)
  if (is.null(base_env)) {
    base_env <- parent.frame()
  }

  helper_env <- new.env(parent = base_env)
  helper_env$ns <- splines::ns
  helper_env$bs <- splines::bs
  environment(formula) <- helper_env
  formula
}

with_formula_helpers_list <- function(formulas) {
  if (is.null(formulas)) {
    return(NULL)
  }
  lapply(formulas, with_formula_helpers)
}
make_tvcqgcomp_config <- function(id,
                                  time_in = NULL,
                                  time_name = NULL,
                                  time_out,
                                  outcome,
                                  outcome_type = "survival",
                                  exposures,
                                  exposure_lags = NULL,
                                  exposure_formulas = NULL,
                                  exposure_types = NULL,
                                  time_varying_covariates,
                                  tvc_lags = NULL,
                                  covtypes = NULL,
                                  passive_tvc_vars = NULL,
                                  time_fixed_covariates,
                                  factor_vars = NULL,
                                  outcome_formula,
                                  tvc_formulas,
                                  q = NULL,
                                  q_levels = NULL,
                                  natural_course = "calibration",
                                  auto_history = TRUE,
                                  baselags = FALSE,
                                  meta_target = "HR",
                                  categorical_sim = "stochastic",
                                  exposome = NULL,
                                  exposure_scale = NULL,
                                  bounded_supports = NULL) {
  normalize_meta_target <- function(x) {
    x <- unique(unlist(lapply(x, function(target) {
      target <- match.arg(target, c("HR", "HD", "RR", "RD", "hazard", "cumrisk_final", "both"))
      if (identical(target, "hazard")) return("HR")
      if (identical(target, "cumrisk_final")) return("RR")
      if (identical(target, "both")) return(c("HR", "RR"))
      target
    }), use.names = FALSE))
    x
  }

  normalize_natural_course <- function(x) {
    x <- unique(unlist(lapply(x, function(mode) {
      mode <- match.arg(mode, c("none", "calibration", "observed", "observed_exposome", "modeled_exposome", "modeled"))
      if (identical(mode, "observed")) return("calibration")
      if (identical(mode, "modeled")) return("modeled_exposome")
      mode
    }), use.names = FALSE))
    x
  }

  normalize_named <- function(x, keys, default = NULL, label = "entries") {
    if (is.null(x)) {
      x <- rep(default, length(keys))
      names(x) <- keys
      return(x)
    }
    if (is.null(names(x))) {
      if (length(x) != length(keys)) {
        stop("Unnamed ", label, " vector must have the same length as its variables.")
      }
      names(x) <- keys
      return(x)
    }
    if (!all(keys %in% names(x))) {
      stop("Missing ", label, " for: ", paste(setdiff(keys, names(x)), collapse = ", "))
    }
    x[keys]
  }

  normalize_bounded_supports <- function(x, valid_names) {
    if (is.null(x)) {
      return(NULL)
    }
    if (!is.list(x) || is.null(names(x))) {
      stop("bounded_supports must be a named list, e.g. list(var = 1:5).")
    }
    out <- x[intersect(names(x), valid_names)]
    for (nm in names(out)) {
      spec <- out[[nm]]
      if (is.null(spec)) next
      if (is.numeric(spec)) {
        vals <- sort(unique(as.numeric(spec)))
      } else if (is.list(spec) && !is.null(spec$values) && is.numeric(spec$values)) {
        vals <- sort(unique(as.numeric(spec$values)))
      } else if (is.list(spec) && !is.null(spec$min) && !is.null(spec$max)) {
        by <- spec$by %||% 1
        vals <- seq(from = as.numeric(spec$min), to = as.numeric(spec$max), by = as.numeric(by))
      } else if (is.character(spec) && length(spec) == 1L && tolower(spec) %in% c("observed", "auto")) {
        vals <- "observed"
      } else {
        stop("Invalid bounded_supports entry for ", nm,
             ". Use numeric values, list(values = ...), list(min=, max=, by=), or 'observed'.")
      }
      if (is.character(vals)) {
        out[[nm]] <- vals
      } else {
        vals <- vals[is.finite(vals)]
        if (length(vals) < 2L) {
          stop("bounded_supports for ", nm, " must contain at least two finite support points.")
        }
        out[[nm]] <- vals
      }
    }
    out
  }

  if (is.null(time_in)) {
    time_in <- time_name
  }
  if (is.null(time_in)) {
    stop("Please supply time_in or time_name.")
  }
  if (is.null(exposure_lags)) {
    exposure_lags <- paste0(exposures, "_lag1")
  }
  if (is.null(tvc_lags)) {
    tvc_lags <- paste0(time_varying_covariates, "_lag1")
  }

  if (is.null(q)) {
    q <- q_levels %||% 4L
  }
  if (is.null(exposome)) {
    exposome <- exposure_scale %||% "quantized"
  }
  bounded_supports <- normalize_bounded_supports(
    bounded_supports,
    valid_names = unique(c(exposures, time_varying_covariates))
  )
  covtypes <- normalize_named(covtypes, time_varying_covariates, default = "categorical", label = "covtype entries")
  exposure_types <- normalize_named(exposure_types, exposures, default = "categorical", label = "exposure type entries")

  if (is.null(passive_tvc_vars)) {
    passive_tvc_vars <- character()
  } else {
    passive_tvc_vars <- unique(as.character(passive_tvc_vars))
    if (!all(passive_tvc_vars %in% time_varying_covariates)) {
      stop(
        "passive_tvc_vars must be a subset of time_varying_covariates. Unknown values: ",
        paste(setdiff(passive_tvc_vars, time_varying_covariates), collapse = ", ")
      )
    }
  }
  active_tvc_vars <- setdiff(time_varying_covariates, passive_tvc_vars)

  if (!is.null(exposure_formulas)) {
    if (!all(exposures %in% names(exposure_formulas))) {
      stop("Missing exposure formulas for: ", paste(setdiff(exposures, names(exposure_formulas)), collapse = ", "))
    }
    exposure_formulas <- exposure_formulas[exposures]
  }

  if (is.null(tvc_formulas)) {
    if (length(active_tvc_vars) > 0L) {
      stop("Missing tvc_formulas for active time-varying covariates: ", paste(active_tvc_vars, collapse = ", "))
    }
    tvc_formulas <- list()
  }
  if (!is.null(names(tvc_formulas))) {
    unknown_tvc_formula_names <- setdiff(names(tvc_formulas), time_varying_covariates)
    if (length(unknown_tvc_formula_names) > 0L) {
      stop("Unknown tvc_formulas entries: ", paste(unknown_tvc_formula_names, collapse = ", "))
    }
  }
  if (!all(active_tvc_vars %in% names(tvc_formulas))) {
    stop("Missing tvc_formulas for active time-varying covariates: ", paste(setdiff(active_tvc_vars, names(tvc_formulas)), collapse = ", "))
  }
  full_tvc_formulas <- stats::setNames(vector("list", length(time_varying_covariates)), time_varying_covariates)
  if (length(tvc_formulas)) {
    full_tvc_formulas[names(tvc_formulas)] <- tvc_formulas[names(tvc_formulas)]
  }
  tvc_formulas <- full_tvc_formulas

  outcome_formula <- with_formula_helpers(outcome_formula)
  tvc_formulas <- with_formula_helpers_list(tvc_formulas)
  exposure_formulas <- with_formula_helpers_list(exposure_formulas)

  q <- as.integer(q)
  exposome <- match.arg(exposome, "quantized")
  categorical_sim <- match.arg(categorical_sim, c("stochastic", "class"))
  outcome_type <- match.arg(outcome_type, "survival")

  list(
    id = id,
    time_in = time_in,
    time_name = time_in,
    time_out = time_out,
    outcome = outcome,
    outcome_type = outcome_type,
    exposures = exposures,
    exposure_lags = exposure_lags,
    exposure_formulas = exposure_formulas,
    exposure_types = exposure_types,
    time_varying_covariates = time_varying_covariates,
    tvc = time_varying_covariates,
    tvc_lags = tvc_lags,
    covtypes = covtypes,
    passive_tvc_vars = passive_tvc_vars,
    active_tvc_vars = active_tvc_vars,
    time_fixed_covariates = time_fixed_covariates,
    baseline_covariates = time_fixed_covariates,
    factor_vars = factor_vars,
    outcome_formula = outcome_formula,
    tvc_formulas = tvc_formulas,
    q = q,
    q_levels = q,
    natural_course = normalize_natural_course(natural_course),
    auto_history = isTRUE(auto_history),
    baselags = isTRUE(baselags),
    meta_target = normalize_meta_target(meta_target),
    categorical_sim = categorical_sim,
    exposome = exposome,
    exposure_scale = exposome,
    bounded_supports = bounded_supports
  )
}

make_no2_tvcqgcomp_config <- function(id = "UniqID",
                                      time_in = "TimeInn",
                                      time_out = "TimeOut",
                                      outcome = "status",
                                      exposures = c("NO2_RF", "BC_RF", "PM1_RF"),
                                      exposure_lags = c("NO2_RF_lag1", "BC_RF_lag1", "PM1_RF_lag1"),
                                      exposure_formulas = NULL,
                                      exposure_types = stats::setNames(rep("categorical", 3L), c("NO2_RF", "BC_RF", "PM1_RF")),
                                      time_varying_covariates = c("CSize", "Urban_form", "ethniconcentration", "deprivation", "dependency", "instability"),
                                      time_fixed_covariates = c(
                                        "age", "sex", "visible_minority", "IndigenousIdentity",
                                        "landed_migration", "marsth", "Education", "employment",
                                        "Occupation", "income_inadequacy"
                                      ),
                                      natural_course = "calibration",
                                      exposome = "quantized") {
  make_tvcqgcomp_config(
    id = id,
    time_in = time_in,
    time_out = time_out,
    outcome = outcome,
    exposures = exposures,
    exposure_lags = exposure_lags,
    exposure_formulas = exposure_formulas,
    exposure_types = exposure_types,
    time_varying_covariates = time_varying_covariates,
    tvc_lags = paste0(time_varying_covariates, "_lag1"),
    covtypes = stats::setNames(rep("categorical", length(time_varying_covariates)), time_varying_covariates),
    time_fixed_covariates = time_fixed_covariates,
    factor_vars = unique(c(
      time_varying_covariates,
      time_fixed_covariates,
      exposures,
      "sex", "visible_minority", "IndigenousIdentity", "landed_migration",
      "marsth", "Education", "employment", "Occupation", "income_inadequacy"
    )),
    outcome_formula = status ~ NO2_RF + BC_RF + PM1_RF +
      TimeOut + I(TimeOut^2) +
      age + sex + visible_minority + IndigenousIdentity + landed_migration +
      marsth + Education + employment + Occupation + income_inadequacy +
      CSize + Urban_form + ethniconcentration + deprivation + dependency + instability,
    tvc_formulas = list(
      CSize = CSize ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        CSize_lag1 + deprivation_lag1 + dependency_lag1 + instability_lag1 + ethniconcentration_lag1,
      Urban_form = Urban_form ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        Urban_form_lag1 + CSize_lag1 + CSize +
        deprivation_lag1 + dependency_lag1 + instability_lag1 + ethniconcentration_lag1,
      ethniconcentration = ethniconcentration ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        ethniconcentration_lag1 + CSize + Urban_form +
        deprivation_lag1 + dependency_lag1 + instability_lag1,
      deprivation = deprivation ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        deprivation_lag1 + CSize + Urban_form + ethniconcentration +
        dependency_lag1 + instability_lag1,
      dependency = dependency ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        dependency_lag1 + CSize + Urban_form + ethniconcentration + deprivation +
        instability_lag1,
      instability = instability ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        instability_lag1 + CSize + Urban_form + ethniconcentration + deprivation + dependency
    ),
    q = 4L,
    natural_course = natural_course,
    auto_history = TRUE,
    meta_target = "HR",
    categorical_sim = "stochastic",
    exposome = exposome
  )
}

quantize_exposures <- function(data, exposures, q = 4L, breaks = NULL) {
  qgcomp::quantize(data = data, expnms = exposures, q = q, breaks = breaks)
}

generate_intervention_scenarios <- function(config, levels = NULL, prefix = "Q") {
  if (is.null(levels)) {
    levels <- 0:(config$q - 1L)
  }
  out <- lapply(levels, function(level) {
    vals <- rep(level, length(config$exposures))
    names(vals) <- config$exposures
    vals
  })
  names(out) <- paste0(prefix, seq_along(out) - 1L)
  out
}

make_intervention_lookup <- function(data,
                                     config,
                                     rules,
                                     id_col = config$id,
                                     time_col = config$time_in) {
  dt <- data.table::as.data.table(data.table::copy(data))
  exposures <- config$exposures

  required_cols <- c(id_col, time_col, exposures)
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols)) {
    stop("Observed data are missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (is.null(rules)) {
    stop("'rules' must be supplied.")
  }
  if (is.atomic(rules) && !is.null(names(rules))) {
    rules <- as.list(rules)
  }
  if (!is.list(rules)) {
    stop("'rules' must be a named list, or a named atomic vector coercible to a list.")
  }
  if (!all(exposures %in% names(rules))) {
    stop("Rules are missing exposures: ", paste(setdiff(exposures, names(rules)), collapse = ", "))
  }

  out <- dt[, c(id_col, time_col), with = FALSE]
  data.table::setnames(out, c(id_col, time_col), c("original_id", time_col))

  n <- nrow(dt)
  for (exp_name in exposures) {
    rule <- rules[[exp_name]]
    values <- if (is.function(rule)) {
      rule(dt[[exp_name]], dt, exp_name)
    } else if (length(rule) == 1L) {
      rep(as.numeric(rule), n)
    } else if (length(rule) == n) {
      as.numeric(rule)
    } else {
      stop(
        "Rule for exposure '", exp_name,
        "' must be a function, a scalar, or a vector of length nrow(data)."
      )
    }
    if (length(values) != n) {
      stop("Rule for exposure '", exp_name, "' returned ", length(values), " values; expected ", n, ".")
    }
    if (any(!is.finite(values))) {
      stop("Rule for exposure '", exp_name, "' returned non-finite values.")
    }
    out[, (exp_name) := values]
  }

  out[]
}

make_q_scenarios <- function(config, levels = NULL, prefix = "Q") {
  generate_intervention_scenarios(config, levels = levels, prefix = prefix)
}

tvcQGComp <- function(...) {
  run_tvcqgcomp(...)
}

run_tvcqgcomp_survival <- function(...) {
  run_tvcqgcomp(...)
}

tvcQGComp_boot <- function(...) {
  bootstrap_tvcqgcomp(...)
}

bootstrap_tvcqgcomp_survival <- function(...) {
  bootstrap_tvcqgcomp(...)
}

tvcQGComp_survival <- function(data,
                               config,
                               mc_size = 10000L,
                               intervention_scenarios = NULL,
                               natural_course = config$natural_course,
                               seed = 1234L,
                               replace_mc = NULL,
                               meta_formula = NULL,
                               meta_family = NULL,
                               meta_target = config$meta_target,
                               fit_meta = NULL,
                               fit_exposure_models = NULL,
                               exposure_scale = config$exposure_scale,
                               quantization_breaks = NULL,
                               joint_effect_fn = NULL,
                               verbose = FALSE) {
  run_tvcqgcomp(
    data = data,
    config = config,
    mc_size = mc_size,
    intervention_scenarios = intervention_scenarios,
    natural_course = natural_course,
    seed = seed,
    replace_mc = replace_mc,
    meta_formula = meta_formula,
    meta_family = meta_family,
    meta_target = meta_target,
    fit_meta = fit_meta,
    fit_exposure_models = fit_exposure_models,
    exposure_scale = exposure_scale,
    quantization_breaks = quantization_breaks,
    joint_effect_fn = joint_effect_fn,
    verbose = verbose
  )
}

tvcQGComp_survival_boot <- function(data,
                                    config,
                                    n_boot = 200L,
                                    mc_size = 10000L,
                                    intervention_scenarios = NULL,
                                    natural_course = config$natural_course,
                                    seed = 1234L,
                                    replace_mc = NULL,
                                    meta_formula = NULL,
                                    meta_family = NULL,
                                    meta_target = config$meta_target,
                                    fit_meta = NULL,
                                    fit_exposure_models = NULL,
                                    exposure_scale = config$exposure_scale,
                                    quantization_breaks = NULL,
                                    joint_effect_fn = NULL,
                                    parallel = FALSE,
                                    n_workers = 2L,
                                    batch_size = NULL,
                                    checkpoint_file = NULL,
                                    resume_from_checkpoint = TRUE,
                                    verbose = TRUE,
                                    stop_on_error = TRUE,
                                    audit_mode = FALSE,
                                    keep_boot_objects = FALSE,
                                    keep_boot_curves = FALSE) {
  bootstrap_tvcqgcomp(
    data = data,
    config = config,
    n_boot = n_boot,
    mc_size = mc_size,
    intervention_scenarios = intervention_scenarios,
    natural_course = natural_course,
    seed = seed,
    replace_mc = replace_mc,
    meta_formula = meta_formula,
    meta_family = meta_family,
    meta_target = meta_target,
    fit_meta = fit_meta,
    fit_exposure_models = fit_exposure_models,
    exposure_scale = exposure_scale,
    quantization_breaks = quantization_breaks,
    joint_effect_fn = joint_effect_fn,
    parallel = parallel,
    n_workers = n_workers,
    batch_size = batch_size,
    checkpoint_file = checkpoint_file,
    resume_from_checkpoint = resume_from_checkpoint,
    verbose = verbose,
    stop_on_error = stop_on_error,
    audit_mode = audit_mode,
    keep_boot_objects = keep_boot_objects,
    keep_boot_curves = keep_boot_curves
  )
}

#' @export
print.tvcqgcomp_run <- function(x,
                                prediction = c("none", "trajectory"),
                                include_natural = FALSE,
                                ...) {
  prediction <- match.arg(prediction)
  cfg <- x$config %||% list()
  natural_modes <- if (!is.null(cfg$natural_course)) paste(cfg$natural_course, collapse = ", ") else "unknown"
  n_scenarios <- if (!is.null(x$interventions)) length(x$interventions) else NA_integer_
  mc_n <- if (!is.null(x$mc_data)) nrow(x$mc_data) else NA_integer_

  cat("tvcQGcomp run\n")
  cat("  Outcome type   :", cfg$outcome_type %||% "survival", "\n")
  cat("  Meta target    :", paste(cfg$meta_target %||% x$meta_target %||% "unknown", collapse = ", "), "\n")
  cat("  Categorical sim:", cfg$categorical_sim %||% "stochastic", "\n")
  cat("  Natural course :", natural_modes, "\n")
  cat("  Scenarios      :", n_scenarios, "\n")
  cat("  MC rows        :", mc_n, "\n")

  if (identical(prediction, "trajectory")) {
    curve_summary <- as_risk_curve_df(
      x,
      include_natural = include_natural,
      natural_label = "Natural Course"
    )
    time_col <- if (!is.null(cfg$time_out) && cfg$time_out %in% names(curve_summary)) {
      cfg$time_out
    } else {
      cfg$time_in %||% "TimeInn"
    }
    pointwise_display <- data.frame(
      "Scenario" = curve_summary$scenario,
      "Time" = curve_summary[[time_col]],
      "Mean Risk" = curve_summary$mean_risk,
      "Mean Survival" = curve_summary$mean_survival,
      check.names = FALSE
    )
    cat("\nRisk trajectory:\n")
    print(pointwise_display, row.names = FALSE)
  }
  if (!is.null(x$meta_effect_summary)) {
    cat("\nMeta-model effect summary:\n")
    effect_display <- data.frame(
      "Meta Target" = x$meta_effect_summary$meta_target,
      "Estimate" = x$meta_effect_summary$estimate,
      "Increment" = x$meta_effect_summary$increment,
      check.names = FALSE
    )
    print(effect_display, row.names = FALSE)
  }

  invisible(x)
}


