`%||%` <- function(x, y) if (is.null(x)) y else x

draw_from_probs_fast <- function(pmat, fallback = "uniform") {
  pmat <- as.matrix(pmat)
  if (is.null(dim(pmat))) {
    pmat <- cbind(pmat)
  }

  pmat[!is.finite(pmat)] <- 0
  pmat <- pmax(pmat, 0)

  rs <- rowSums(pmat)
  bad <- rs <= 0
  if (any(bad)) {
    if (identical(fallback, "uniform")) {
      pmat[bad, ] <- 1
      rs[bad] <- ncol(pmat)
    } else {
      stop("Bad probability rows in draw_from_probs_fast().")
    }
  }

  pmat <- pmat / rs
  cs <- pmat
  if (ncol(pmat) > 1L) {
    for (j in 2:ncol(pmat)) {
      cs[, j] <- cs[, j] + cs[, j - 1L]
    }
  }
  cs[, ncol(pmat)] <- 1
  u <- stats::runif(nrow(cs))
  max.col(u <= cs, ties.method = "first")
}

safe_predict_probs <- function(model, newdata, varname, time_value, id_value = NULL) {
  p <- tryCatch(
    predict_tvcqgcomp_model(model, newdata = newdata, type = "probs"),
    error = function(e) {
      stop(
        "Predict error for ", varname,
        " at time ", time_value,
        if (!is.null(id_value)) paste0(" (id=", id_value, ")") else "",
        ": ", conditionMessage(e)
      )
    }
  )

  p <- as.matrix(p)
  bad_rows <- !is.finite(rowSums(p, na.rm = TRUE)) | rowSums(p, na.rm = TRUE) <= 0
  if (any(bad_rows)) {
    stop("Invalid probability rows for ", varname, " at time ", time_value)
  }
  p
}

safe_predict_response <- function(model, newdata, varname, time_value, id_value = NULL) {
  pred <- tryCatch(
    as.numeric(predict_tvcqgcomp_model(model, newdata = newdata, type = "response")),
    error = function(e) {
      stop(
        "Predict error for ", varname,
        " at time ", time_value,
        if (!is.null(id_value)) paste0(" (id=", id_value, ")") else "",
        ": ", conditionMessage(e)
      )
    }
  )
  bad <- which(!is.finite(pred) | is.na(pred))
  if (length(bad)) {
    bad_newdata <- data.table::as.data.table(newdata[bad, , drop = FALSE])
    missing_cols <- names(bad_newdata)[vapply(
      bad_newdata,
      function(x) any(!is.finite(suppressWarnings(as.numeric(x))) | is.na(x)),
      logical(1)
    )]
    missing_cols <- setdiff(
      missing_cols,
      names(bad_newdata)[vapply(bad_newdata, is.character, logical(1))]
    )
    missing_msg <- if (length(missing_cols)) {
      paste0(" Missing/invalid predictors in: ", paste(utils::head(missing_cols, 10L), collapse = ", "),
             if (length(missing_cols) > 10L) ", ..." else "")
    } else {
      ""
    }
    stop(
      "Predict returned non-finite values for ", varname,
      " at time ", time_value,
      if (!is.null(id_value)) paste0(" (id=", id_value, ")") else "",
      ". Problem rows: ", length(bad), ".", missing_msg
    )
  }
  pred
}

assert_simulated_values <- function(values, value_name, time_value, id_value = NULL) {
  bad <- which(is.na(values) | !is.finite(values))
  if (!length(bad)) {
    return(invisible(values))
  }
  stop(
    "Simulation produced non-finite values for ", value_name,
    " at time ", time_value,
    if (!is.null(id_value)) paste0(" (id=", id_value, ")") else "",
    ". Problem rows: ", length(bad), "."
  )
}

fit_gaussian_support <- function(model, observed) {
  list(
    sigma = stats::sd(stats::residuals(model, type = "response"), na.rm = TRUE),
    min = min(observed, na.rm = TRUE),
    max = max(observed, na.rm = TRUE)
  )
}

fit_bounded_normal_model <- function(formula, data, value_name) {
  observed_raw <- data[[value_name]]
  observed <- suppressWarnings(as.numeric(observed_raw))
  if (any(!is.na(observed_raw) & is.na(observed))) {
    stop("Bounded normal variable ", value_name, " must be numeric/integer (or factor coercible to numeric).")
  }
  min_val <- min(observed, na.rm = TRUE)
  max_val <- max(observed, na.rm = TRUE)
  if (!is.finite(min_val) || !is.finite(max_val) || max_val <= min_val) {
    stop("Bounded normal variable ", value_name, " must have finite support with max > min.")
  }

  normalized <- (observed - min_val) / (max_val - min_val)
  fit_data <- data.table::copy(data)
  fit_data[[value_name]] <- normalized
  model <- stats::glm(formula, data = fit_data, family = stats::gaussian())
  sigma <- stats::sd(stats::residuals(model, type = "response"), na.rm = TRUE)
  if (!is.finite(sigma) || is.na(sigma)) sigma <- 0

  structure(
    list(
      model = model,
      min = min_val,
      max = max_val,
      sigma = sigma,
      value_name = value_name,
      snap = FALSE,
      support_values = NULL,
      original_is_integer = is.integer(observed_raw)
    ),
    class = "tvcqgcomp_bounded_normal"
  )
}

resolve_support_values <- function(observed, support_spec = NULL) {
  if (is.null(support_spec) || (is.character(support_spec) && length(support_spec) == 1L &&
                                tolower(support_spec) %in% c("observed", "auto"))) {
    vals <- sort(unique(observed[is.finite(observed)]))
  } else {
    vals <- as.numeric(support_spec)
    vals <- vals[is.finite(vals)]
    vals <- sort(unique(vals))
  }
  if (length(vals) < 2L) {
    stop("Support for snapped bounded normal variables must contain at least two finite values.")
  }
  vals
}

snap_to_support <- function(values, support_values) {
  values_num <- as.numeric(unlist(values, use.names = FALSE))
  support_num <- as.numeric(unlist(support_values, use.names = FALSE))
  support_num <- sort(unique(support_num[is.finite(support_num)]))
  if (!length(support_num)) {
    stop("snap_to_support() received an empty or non-finite support grid.")
  }
  if (length(support_num) == 1L) {
    return(rep(support_num[[1L]], length(values_num)))
  }
  idx <- vapply(
    values_num,
    function(z) {
      if (!is.finite(z) || is.na(z)) {
        return(NA_integer_)
      }
      which.min(abs(z - support_num))
    },
    integer(1)
  )
  out <- rep(NA_real_, length(values_num))
  ok <- !is.na(idx)
  out[ok] <- support_num[idx[ok]]
  out
}

simulate_binary <- function(model, newdata) {
  p <- safe_predict_response(model, newdata, varname = "binary", time_value = "unknown")
  stats::rbinom(length(p), size = 1L, prob = pmin(pmax(p, 0), 1))
}

simulate_continuous <- function(model, newdata, support, value_name = NULL, time_value = NULL, id_value = NULL) {
  pred <- safe_predict_response(
    model,
    newdata,
    varname = value_name %||% "continuous",
    time_value = time_value %||% "unknown",
    id_value = id_value
  )
  vals <- pred + stats::rnorm(length(pred), sd = support$sigma)
  pmax(pmin(vals, support$max), support$min)
}

simulate_bounded_normal <- function(model, newdata, value_name = NULL, time_value = NULL, id_value = NULL) {
  pred <- safe_predict_response(
    model,
    newdata,
    varname = value_name %||% model$value_name %||% "bounded_normal",
    time_value = time_value %||% "unknown",
    id_value = id_value
  )
  normalized <- pred + stats::rnorm(length(pred), sd = model$sigma)
  vals <- normalized * (model$max - model$min) + model$min
  vals <- pmax(pmin(vals, model$max), model$min)
  if (isTRUE(model$snap)) {
    vals <- snap_to_support(vals, model$support_values)
    if (isTRUE(model$original_is_integer)) {
      vals <- as.integer(round(vals))
    }
  }
  vals
}

force_levels_dt <- function(dt, reference_dt, vars) {
  for (v in vars) {
    if (v %in% names(dt) && v %in% names(reference_dt) && is.factor(reference_dt[[v]])) {
      dt[, (v) := factor(as.character(get(v)), levels = levels(reference_dt[[v]]), ordered = is.ordered(reference_dt[[v]]))]
    }
  }
  dt[]
}

standardize_output_columns <- function(dt, config, extra_cols = NULL) {
  id_cols <- intersect(c("uid", "sim_person", "original_id", config$id, config$time_in, config$time_out), names(dt))
  base_cols <- intersect(config$baseline_covariates, names(dt))
  tvc_cols <- intersect(config$tvc, names(dt))
  tvc_lag_cols <- intersect(config$tvc_lags, names(dt))
  exposure_lag_cols <- intersect(config$exposure_lags, names(dt))
  exposure_cols <- intersect(config$exposures, names(dt))
  outcome_cols <- intersect(config$outcome, names(dt))
  extra_cols <- intersect(extra_cols, names(dt))

  ordered_cols <- unique(c(id_cols, base_cols, tvc_cols, tvc_lag_cols, exposure_lag_cols, exposure_cols, outcome_cols, extra_cols))
  data.table::setcolorder(dt, c(ordered_cols, setdiff(names(dt), ordered_cols)))
  dt[]
}
get_required_simulation_columns <- function(config, extra_cols = NULL) {
  tvc_formulas <- Filter(Negate(is.null), config$tvc_formulas)
  exposure_formulas <- Filter(Negate(is.null), config$exposure_formulas %||% list())
  formula_vars <- unique(c(
    all.vars(config$outcome_formula),
    if (length(tvc_formulas)) unlist(lapply(tvc_formulas, all.vars), use.names = FALSE) else character(),
    if (length(exposure_formulas)) unlist(lapply(exposure_formulas, all.vars), use.names = FALSE) else character()
  ))

  unique(c(
    "uid", "sim_person", "original_id",
    config$id,
    config$time_in,
    config$time_out,
    config$outcome,
    config$baseline_covariates,
    config$tvc,
    config$tvc_lags,
    config$exposures,
    config$exposure_lags,
    formula_vars,
    extra_cols
  ))
}

trim_simulation_columns <- function(dt, config, extra_cols = NULL) {
  keep_cols <- intersect(get_required_simulation_columns(config, extra_cols = extra_cols), names(dt))
  dt <- dt[, ..keep_cols]
  standardize_output_columns(dt, config, extra_cols = extra_cols)
}
prepare_exposure_scale_data <- function(data, config, exposure_scale = config$exposure_scale, q = config$q_levels, breaks = NULL) {
  dt <- data.table::as.data.table(data.table::copy(data))
  exposure_scale <- match.arg(exposure_scale, "quantized")

  qprep <- quantize_exposures(data = as.data.frame(dt), exposures = config$exposures, q = q, breaks = breaks)
  qdt <- data.table::as.data.table(qprep$data)

  # Existing lag columns remain on the raw exposure scale after quantization.
  # Drop them so prepare_tvcqgcomp_data() rebuilds lag histories from the
  # quantized exposures before fitting any models.
  lag_cols <- intersect(config$exposure_lags, names(qdt))
  if (length(lag_cols) > 0L) {
    qdt[, (lag_cols) := NULL]
  }

  list(data = qdt, breaks = qprep$breaks)
}

fit_single_stochastic_model <- function(formula, data, model_type, value_name = NULL, support_spec = NULL) {
  model_type <- tolower(model_type)
  if (model_type %in% c("categorical", "multinomial")) {
    fit <- nnet::multinom(formula, data = data, trace = FALSE, Hess = FALSE, MaxNWts = 5000)
  } else if (model_type %in% c("ordinal", "ordinal_logit", "ordinal_polr", "ordered_logit", "polr")) {
    if (is.null(value_name)) stop("value_name is required for ordinal models.")
    fit_data <- data.table::copy(data)
    y <- fit_data[[value_name]]
    if (is.ordered(y)) {
      y_ord <- y
    } else if (is.factor(y)) {
      y_ord <- factor(y, levels = levels(y), ordered = TRUE)
    } else {
      y_chr <- as.character(y)
      y_num <- suppressWarnings(as.numeric(y_chr))
      if (!all(is.na(y_chr) | is.finite(y_num))) {
        lev <- sort(unique(stats::na.omit(y_chr)))
      } else {
        lev <- as.character(sort(unique(stats::na.omit(y_num))))
      }
      y_ord <- factor(y_chr, levels = lev, ordered = TRUE)
    }
    fit_data[[value_name]] <- y_ord
    fit <- MASS::polr(formula, data = fit_data, method = "logistic", Hess = FALSE, model = FALSE)
  } else if (model_type %in% c("binary", "binomial")) {
    fit <- stats::glm(formula, data = data, family = stats::binomial())
  } else if (model_type %in% c("normal", "continuous")) {
    fit <- stats::lm(formula, data = data)
  } else if (model_type %in% c("bounded normal", "bounded_normal", "bounded normal snap", "bounded_normal_snap", "bounded normal ordinal", "bounded_normal_ordinal")) {
    if (is.null(value_name)) stop("value_name is required for bounded normal models.")
    fit <- fit_bounded_normal_model(formula, data, value_name = value_name)
    if (model_type %in% c("bounded normal snap", "bounded_normal_snap", "bounded normal ordinal", "bounded_normal_ordinal")) {
      observed <- suppressWarnings(as.numeric(data[[value_name]]))
      fit$snap <- TRUE
      fit$support_values <- resolve_support_values(observed, support_spec = support_spec)
      fit$min <- min(fit$support_values, na.rm = TRUE)
      fit$max <- max(fit$support_values, na.rm = TRUE)
    }
  } else {
    stop("Unsupported model type: ", model_type)
  }
  attach_prediction_support(fit, formula, data)
}

simulate_from_model <- function(model, newdata, model_type, support = NULL, value_name = NULL, time_value = NULL, id_value = NULL, categorical_sim = "stochastic") {
  model_type <- tolower(model_type)
  if (model_type %in% c("categorical", "multinomial", "ordinal", "ordinal_logit", "ordinal_polr", "ordered_logit", "polr")) {
    categorical_sim <- match.arg(categorical_sim, c("stochastic", "class"))
    if (identical(categorical_sim, "class")) {
      pred_class <- tryCatch(
        predict_tvcqgcomp_model(model, newdata = newdata, type = "class"),
        error = function(e) {
          stop(
            "Predict error for ", value_name %||% "value",
            " at time ", time_value %||% "unknown",
            if (!is.null(id_value)) paste0(" (id=", id_value, ")") else "",
            ": ", conditionMessage(e)
          )
        }
      )
      if (model_type %in% c("ordinal", "ordinal_logit", "ordinal_polr", "ordered_logit", "polr") &&
          value_name %in% names(newdata) && (is.numeric(newdata[[value_name]]) || is.integer(newdata[[value_name]]))) {
        return(suppressWarnings(as.numeric(as.character(pred_class))))
      }
      if (value_name %in% names(newdata) && is.factor(newdata[[value_name]])) {
        return(factor(as.character(pred_class), levels = levels(newdata[[value_name]]), ordered = is.ordered(newdata[[value_name]])))
      }
      return(pred_class)
    }
    p <- safe_predict_probs(model, newdata, value_name %||% "value", time_value %||% "unknown", id_value = id_value)
    p <- as.matrix(p)
    level_names <- colnames(p)
    if (!is.null(level_names) && value_name %in% names(newdata) && is.factor(newdata[[value_name]])) {
      miss <- setdiff(levels(newdata[[value_name]]), level_names)
      for (m in miss) {
        p <- cbind(p, stats::setNames(rep(0, nrow(p)), m))
      }
      p <- p[, levels(newdata[[value_name]]), drop = FALSE]
      k <- draw_from_probs_fast(p)
      return(factor(levels(newdata[[value_name]])[k], levels = levels(newdata[[value_name]]), ordered = is.ordered(newdata[[value_name]])))
    }
    k <- draw_from_probs_fast(p)
    if (!is.null(level_names)) {
      vals <- level_names[k]
      if (model_type %in% c("ordinal", "ordinal_logit", "ordinal_polr", "ordered_logit", "polr")) {
        if (value_name %in% names(newdata) && (is.numeric(newdata[[value_name]]) || is.integer(newdata[[value_name]]))) {
          return(suppressWarnings(as.numeric(vals)))
        }
        return(factor(vals, levels = level_names, ordered = TRUE))
      }
      if (value_name %in% names(newdata) && (is.numeric(newdata[[value_name]]) || is.integer(newdata[[value_name]]))) {
        return(suppressWarnings(as.numeric(vals)))
      }
      return(vals)
    }
    return(k - 1L)
  }
  if (model_type %in% c("binary", "binomial")) {
    pred <- safe_predict_response(
      model,
      newdata,
      varname = value_name %||% "binary",
      time_value = time_value %||% "unknown",
      id_value = id_value
    )
    return(stats::rbinom(length(pred), size = 1L, prob = pmin(pmax(pred, 0), 1)))
  }
  if (model_type %in% c("normal", "continuous")) {
    return(simulate_continuous(
      model,
      newdata,
      support,
      value_name = value_name,
      time_value = time_value,
      id_value = id_value
    ))
  }
  if (model_type %in% c("bounded normal", "bounded_normal", "bounded normal snap", "bounded_normal_snap", "bounded normal ordinal", "bounded_normal_ordinal")) {
    return(simulate_bounded_normal(
      model,
      newdata,
      value_name = value_name,
      time_value = time_value,
      id_value = id_value
    ))
  }
  stop("Unsupported model type: ", model_type)
}

fit_tvcqgcomp_models <- function(data, config, progress_fn = NULL, fit_exposure_models = TRUE) {
  dt <- prepare_tvcqgcomp_data(data, config)
  outcome_dt <- dt
  outcome_formula <- config$outcome_formula

  outcome_model <- speedglm::speedglm(
    formula = outcome_formula,
    data = outcome_dt,
    family = stats::binomial("logit"),
    sparse = FALSE
  )
  if (is.function(progress_fn)) progress_fn("fit outcome model")

  tvc_models <- stats::setNames(vector("list", length(config$tvc)), config$tvc)
  tvc_support <- stats::setNames(vector("list", length(config$tvc)), config$tvc)
  for (v in config$tvc) {
    if (v %in% (config$passive_tvc_vars %||% character())) {
      tvc_models[[v]] <- NULL
      tvc_support[[v]] <- NULL
      if (is.function(progress_fn)) progress_fn(paste0("skip passive confounder model: ", v))
      next
    }
    model <- fit_single_stochastic_model(
      config$tvc_formulas[[v]],
      dt,
      config$covtypes[[v]],
      value_name = v,
      support_spec = config$bounded_supports[[v]] %||% NULL
    )
    tvc_models[[v]] <- model
    tvc_support[[v]] <- if (tolower(config$covtypes[[v]]) %in% c("normal", "continuous")) fit_gaussian_support(model, dt[[v]]) else NULL
    if (is.function(progress_fn)) progress_fn(paste0("fit confounder model: ", v))
  }

  exposure_models <- NULL
  exposure_support <- NULL
  if (isTRUE(fit_exposure_models) && !is.null(config$exposure_formulas)) {
    exposure_models <- stats::setNames(vector("list", length(config$exposures)), config$exposures)
    exposure_support <- stats::setNames(vector("list", length(config$exposures)), config$exposures)
    for (v in config$exposures) {
      model <- fit_single_stochastic_model(
        config$exposure_formulas[[v]],
        dt,
        config$exposure_types[[v]],
        value_name = v,
        support_spec = config$bounded_supports[[v]] %||% NULL
      )
      exposure_models[[v]] <- model
      exposure_support[[v]] <- if (tolower(config$exposure_types[[v]]) %in% c("normal", "continuous")) fit_gaussian_support(model, dt[[v]]) else NULL
      if (is.function(progress_fn)) progress_fn(paste0("fit exposure model: ", v))
    }
  }

  list(
    outcome_model = outcome_model,
    outcome_formula = outcome_formula,
    tvc_models = tvc_models,
    tvc_support = tvc_support,
    exposure_models = exposure_models,
    exposure_support = exposure_support
  )
}
create_mc_skeleton <- function(data, config, mc_size = 10000L, seed = 1234L, replace = NULL) {
  dt <- prepare_tvcqgcomp_data(data, config)
  force_levels_dt(dt, dt, intersect(config$factor_vars, names(dt)))
  history_plan <- get_history_plan(config)
  history_lags <- if (nrow(history_plan$lag_specs)) history_plan$lag_specs$lag_name else character()

  set.seed(seed)

  base_time <- min(dt[[config$time_in]], na.rm = TRUE)
  times <- sort(unique(dt[[config$time_in]]))
  n_times <- length(times)

  tvc_levels <- lapply(config$tvc, function(v) if (is.factor(dt[[v]])) levels(dt[[v]]) else NULL)
  names(tvc_levels) <- config$tvc

  base_cols <- unique(c(config$id, config$baseline_covariates))
  if ("original_id" %in% names(dt)) {
    base_cols <- unique(c(base_cols, "original_id"))
  }
  base_dt <- dt[, ..base_cols][, .SD[1L], by = c(config$id)]
  base_ids <- base_dt[[config$id]]

  if (is.null(replace)) {
    replace <- mc_size > length(base_ids)
  }

  t0_cols <- unique(c(config$id, config$tvc))
  base_tvc <- dt[get(config$time_in) == base_time, ..t0_cols]
  data.table::setnames(base_tvc, config$tvc, paste0(config$tvc, '_t0'))

  data.table::setkeyv(base_dt, config$id)
  data.table::setkeyv(base_tvc, config$id)
  base_dt <- base_tvc[base_dt]

  draw_ids <- sample(base_ids, size = mc_size, replace = replace)
  sim_base <- base_dt[J(draw_ids), nomatch = 0L, allow.cartesian = TRUE]
  sim_base[, sim_person := .I]
  sim_base[, uid := paste0('sim_', sim_person)]
  if (!("original_id" %in% names(sim_base))) {
    sim_base[, original_id := get(config$id)]
  }

  sim_dt <- sim_base[rep(seq_len(.N), each = n_times)]
  sim_dt[, (config$time_in) := rep(times, times = mc_size)]

  time_lookup <- unique(dt[, c(config$time_in, config$time_out), with = FALSE])
  data.table::setnames(time_lookup, c(config$time_in, config$time_out), c('.time_in', '.time_out'))
  sim_dt <- merge(sim_dt, time_lookup, by.x = config$time_in, by.y = '.time_in', all.x = TRUE, sort = FALSE)
  data.table::setnames(sim_dt, '.time_out', config$time_out)

  exposure_cols <- unique(c(config$id, config$time_in, config$exposures, intersect(history_lags, names(dt))))
  exposure_dt <- dt[, ..exposure_cols]
  data.table::setkeyv(exposure_dt, c(config$id, config$time_in))

  passive_tvc_vars <- config$passive_tvc_vars %||% character()
  passive_cols <- unique(c(config$id, config$time_in, passive_tvc_vars))
  passive_dt <- NULL
  if (length(passive_tvc_vars) > 0L) {
    passive_dt <- dt[, ..passive_cols]
    data.table::setkeyv(passive_dt, c(config$id, config$time_in))
  }

  data.table::setkeyv(sim_dt, c(config$id, config$time_in))
  sim_dt <- exposure_dt[sim_dt]
  if (!is.null(passive_dt)) {
    sim_dt <- passive_dt[sim_dt]
  }
  sim_dt[, (config$id) := NULL]



  for (v in config$tvc) {
    if (v %in% passive_tvc_vars) {
      next
    }
    vtype <- tolower(config$covtypes[[v]])
    if (vtype %in% c('categorical', 'multinomial')) {
      sim_dt[, (v) := factor(NA, levels = tvc_levels[[v]])]
    } else if (vtype %in% c('ordinal', 'ordinal_logit', 'ordinal_polr', 'ordered_logit', 'polr')) {
      if (is.factor(dt[[v]])) {
        sim_dt[, (v) := factor(NA, levels = levels(dt[[v]]), ordered = TRUE)]
      } else if (is.integer(dt[[v]])) {
        sim_dt[, (v) := NA_integer_]
      } else {
        sim_dt[, (v) := NA_real_]
      }
    } else if (vtype %in% c('binary', 'binomial')) {
      sim_dt[, (v) := NA_integer_]
    } else {
      sim_dt[, (v) := NA_real_]
    }
  }

  for (v in config$tvc) {
    if (v %in% passive_tvc_vars) {
      next
    }
    sim_dt[get(config$time_in) == base_time, (v) := get(paste0(v, '_t0'))]
  }

  t0_names <- intersect(paste0(config$tvc, '_t0'), names(sim_dt))
  if (length(t0_names) > 0L) {
    for (nm in t0_names) {
      sim_dt[, (nm) := NULL]
    }
  }

  if (nrow(history_plan$lag_specs)) {
    for (i in seq_len(nrow(history_plan$lag_specs))) {
      lag_name <- history_plan$lag_specs$lag_name[[i]]
      base_name <- history_plan$lag_specs$base_name[[i]]
      if (lag_name %in% names(sim_dt) || !(base_name %in% names(dt))) next

      if (is.factor(dt[[base_name]])) {
        sim_dt[, (lag_name) := factor(NA, levels = levels(dt[[base_name]]), ordered = is.ordered(dt[[base_name]]))]
      } else if (is.integer(dt[[base_name]])) {
        sim_dt[, (lag_name) := NA_integer_]
      } else {
        sim_dt[, (lag_name) := NA_real_]
      }
    }
  }

  sim_dt[, (config$outcome) := NA_real_]
  data.table::setorderv(sim_dt, c('uid', config$time_in))
  trim_simulation_columns(sim_dt, config)
}

populate_history_for_time <- function(dt, config, history_plan, time_value, times) {
  if (!nrow(history_plan$lag_specs)) {
    return(dt)
  }

  pos <- match(time_value, times)
  if (is.na(pos)) {
    return(dt)
  }

  for (i in seq_len(nrow(history_plan$lag_specs))) {
    lag_name <- history_plan$lag_specs$lag_name[[i]]
    base_name <- history_plan$lag_specs$base_name[[i]]
    lag_n <- history_plan$lag_specs$lag_n[[i]]

    if (!(lag_name %in% names(dt)) || !(base_name %in% names(dt))) next
    if (pos <= lag_n) {
      if (isTRUE(config$baselags)) {
        ref_time <- times[[1L]]
        ref <- dt[get(config$time_in) == ref_time, .(uid, .history_value__ = get(base_name))]
        data.table::setkeyv(ref, "uid")
        dt[get(config$time_in) == time_value, (lag_name) := ref[.SD, on = "uid", .history_value__]]
      } else if (is.factor(dt[[base_name]])) {
        dt[get(config$time_in) == time_value, (lag_name) := factor(
          levels(dt[[base_name]])[[1L]],
          levels = levels(dt[[base_name]]),
          ordered = is.ordered(dt[[base_name]])
        )]
      } else if (is.integer(dt[[base_name]])) {
        dt[get(config$time_in) == time_value, (lag_name) := 0L]
      } else if (is.numeric(dt[[base_name]])) {
        dt[get(config$time_in) == time_value, (lag_name) := 0]
      }
      next
    }

    ref_time <- times[[pos - lag_n]]
    ref <- dt[get(config$time_in) == ref_time, .(uid, .history_value__ = get(base_name))]
    data.table::setkeyv(ref, "uid")
    dt[get(config$time_in) == time_value, (lag_name) := ref[.SD, on = "uid", .history_value__]]
  }

  dt
}

carry_forward_passive_tvcs <- function(dt, config, time_value, prev_tvc = NULL) {
  passive_tvc_vars <- config$passive_tvc_vars %||% character()
  if (!length(passive_tvc_vars)) {
    return(dt)
  }

  idx_expr <- dt[[config$time_in]] == time_value
  if (!any(idx_expr)) {
    return(dt)
  }

  if (is.null(prev_tvc)) {
    prev_time <- max(dt[[config$time_in]][dt[[config$time_in]] < time_value], na.rm = TRUE)
    if (!is.finite(prev_time)) {
      return(dt)
    }
    prev_cols <- c("uid", passive_tvc_vars)
    prev_tvc <- dt[get(config$time_in) == prev_time, ..prev_cols]
    data.table::setnames(prev_tvc, passive_tvc_vars, paste0(passive_tvc_vars, "_prev"))
    data.table::setkeyv(prev_tvc, "uid")
  }

  for (v in passive_tvc_vars) {
    prev_name <- paste0(v, "_prev")
    if (!(prev_name %in% names(prev_tvc)) || !(v %in% names(dt))) next

    fill_vals <- prev_tvc[dt[idx_expr], on = "uid", get(prev_name)]
    cur_vals <- dt[[v]][idx_expr]
    missing_idx <- is.na(cur_vals)
    if (!any(missing_idx)) next

    replace_vals <- fill_vals[missing_idx]
    if (is.factor(dt[[v]])) {
      cur_chr <- as.character(cur_vals)
      cur_chr[missing_idx] <- as.character(replace_vals)
      dt[idx_expr, (v) := factor(cur_chr, levels = levels(dt[[v]]), ordered = is.ordered(dt[[v]]))]
    } else if (is.integer(dt[[v]])) {
      cur_vals[missing_idx] <- as.integer(replace_vals)
      dt[idx_expr, (v) := as.integer(cur_vals)]
    } else if (is.numeric(dt[[v]])) {
      cur_vals[missing_idx] <- as.numeric(replace_vals)
      dt[idx_expr, (v) := as.numeric(cur_vals)]
    } else {
      cur_vals[missing_idx] <- replace_vals
      dt[idx_expr, (v) := cur_vals]
    }
  }

  dt
}

predict_observed_natural_course <- function(data, config, model_fit) {
  dt <- prepare_tvcqgcomp_data(data, config)
  force_levels_dt(dt, dt, intersect(config$factor_vars, names(dt)))
  dt[, hazard := predict_tvcqgcomp_model(model_fit$outcome_model, newdata = dt, type = "response")]
  dt[, hazard := pmin(pmax(hazard, 0), 1)]
  dt[, cumrisk := 1 - cumprod(1 - hazard), by = c(config$id)]
  dt[, survival := 1 - cumrisk]
  trim_simulation_columns(dt, config, extra_cols = c("hazard", "cumrisk", "survival"))
}

predict_natural_course <- function(data,
                                   config,
                                   model_fit,
                                   mc_data = NULL,
                                   mode = c("calibration", "observed", "observed_exposome", "modeled_exposome", "modeled")) {
  mode <- match.arg(mode)
  mode <- normalize_natural_course_modes(mode)[[1L]]
  if (identical(mode, "calibration")) {
    return(predict_observed_natural_course(data, config, model_fit))
  }
  if (identical(mode, "observed_exposome")) {
    if (is.null(mc_data)) {
      stop("mc_data is required for observed_exposome natural course prediction.")
    }
    return(simulate_observed_exposome_natural_course(mc_data, config, model_fit))
  }
  if (is.null(mc_data)) {
    stop("mc_data is required for modeled natural course prediction.")
  }
  if (is.null(model_fit$exposure_models)) {
    stop("Modeled natural course requires exposure_formulas and fitted exposure models in the config.")
  }
  simulate_modeled_exposome_natural_course(mc_data, config, model_fit)
}

normalize_scenario <- function(scenario_values, config, scenario_name = NULL) {
  exposures <- config$exposures
  if (is.function(scenario_values)) {
    return(list(
      type = "rule",
      fun = scenario_values,
      trace_env = NULL,
      trace_name = "rule_trace"
    ))
  }
  if (is.list(scenario_values) && identical(scenario_values$type %||% NULL, "lookup")) {
    return(scenario_values)
  }
  if (is.list(scenario_values) && identical(scenario_values$type %||% NULL, "rule")) {
    fun <- scenario_values$fun %||% scenario_values$rule %||% NULL
    if (!is.function(fun)) {
      stop("Rule scenario ", scenario_name %||% "", " must provide a callable 'fun' or 'rule'.")
    }
    return(list(
      type = "rule",
      fun = fun,
      trace_env = scenario_values$trace_env %||% NULL,
      trace_name = scenario_values$trace_name %||% "rule_trace"
    ))
  }
  if (is.list(scenario_values) && identical(scenario_values$type %||% NULL, "static")) {
    vals <- scenario_values$values %||% NULL
    if (is.null(vals) || !all(exposures %in% names(vals))) {
      stop("Normalized static scenario ", scenario_name %||% "", " is missing exposure values.")
    }
    vals <- as.numeric(vals[exposures])
    names(vals) <- exposures
    return(list(type = "static", values = vals))
  }
  if (data.table::is.data.table(scenario_values) || is.data.frame(scenario_values)) {
    lookup <- data.table::as.data.table(data.table::copy(scenario_values))
    id_candidates <- c("original_id", config$id)
    id_col <- id_candidates[id_candidates %in% names(lookup)][1L]
    if (is.na(id_col) || is.null(id_col)) {
      stop("Lookup scenario ", scenario_name %||% "", " must include either 'original_id' or '", config$id, "'.")
    }
    if (!(config$time_in %in% names(lookup))) {
      stop("Lookup scenario ", scenario_name %||% "", " must include the time column '", config$time_in, "'.")
    }
    if (!all(exposures %in% names(lookup))) {
      stop("Lookup scenario ", scenario_name %||% "", " is missing exposure columns: ", paste(setdiff(exposures, names(lookup)), collapse = ", "))
    }
    lookup <- lookup[, c(id_col, config$time_in, exposures), with = FALSE]
    data.table::setnames(lookup, id_col, "original_id")
    dup_n <- lookup[, .N, by = c("original_id", config$time_in)][N > 1L, .N]
    if (dup_n > 0L) {
      stop("Lookup scenario ", scenario_name %||% "", " has duplicate original_id-time rows.")
    }
    return(list(type = "lookup", data = lookup))
  }
  if (length(scenario_values) == 1L) {
    vals <- rep(as.numeric(scenario_values), length(exposures))
    names(vals) <- exposures
    return(list(type = "static", values = vals))
  }
  if (is.null(names(scenario_values))) {
    if (length(scenario_values) != length(exposures)) {
      stop("Scenario ", scenario_name %||% "", " must have length 1 or match the number of exposures.")
    }
    vals <- as.numeric(scenario_values)
    names(vals) <- exposures
    return(list(type = "static", values = vals))
  }
  if (!all(exposures %in% names(scenario_values))) {
    stop("Scenario ", scenario_name %||% "", " is missing exposure values for: ", paste(setdiff(exposures, names(scenario_values)), collapse = ", "))
  }
  vals <- as.numeric(scenario_values[exposures])
  names(vals) <- exposures
  list(type = "static", values = vals)
}

apply_lookup_scenario_at_time <- function(dt, scenario_spec, config, time_value) {
  idx <- dt[[config$time_in]] == time_value
  if (!any(idx)) {
    return(dt)
  }
  join_dt <- scenario_spec$data[get(config$time_in) == time_value]
  if (!nrow(join_dt)) {
    stop("Lookup scenario is missing rows for time ", time_value, ".")
  }
  cur_dt <- dt[idx]
  matched <- join_dt[cur_dt, on = c("original_id", config$time_in)]
  for (exp_name in config$exposures) {
    vals <- matched[[exp_name]]
    if (any(is.na(vals))) {
      stop("Lookup scenario has missing assigned values for exposure '", exp_name, "' at time ", time_value, ".")
    }
    dt[idx, (exp_name) := vals]
  }
  dt
}

apply_rule_scenario_at_time <- function(dt, scenario_spec, config, time_value) {
  idx <- dt[[config$time_in]] == time_value
  if (!any(idx)) {
    return(dt)
  }
  cur_dt <- data.table::copy(dt[idx])
  before_stats <- cur_dt[, c(
    list(TimeInn = unique(get(config$time_in)), n = .N),
    stats::setNames(
      as.list(unlist(lapply(config$exposures, function(exp_name) mean(get(exp_name), na.rm = TRUE)))),
      paste0(config$exposures, "_before")
    )
  )]
  out_dt <- scenario_spec$fun(cur_dt, config = config, time_value = time_value)
  out_dt <- data.table::as.data.table(out_dt)
  for (exp_name in config$exposures) {
    if (!(exp_name %in% names(out_dt))) {
      stop("Rule scenario did not return exposure column '", exp_name, "' at time ", time_value, ".")
    }
    vals <- out_dt[[exp_name]]
    if (length(vals) != sum(idx)) {
      stop("Rule scenario returned ", length(vals), " values for exposure '", exp_name,
           "' at time ", time_value, " but expected ", sum(idx), ".")
    }
    if (any(!is.finite(vals) | is.na(vals))) {
      stop("Rule scenario produced missing/non-finite values for exposure '", exp_name,
           "' at time ", time_value, ".")
    }
    dt[idx, (exp_name) := vals]
  }
  trace_env <- scenario_spec$trace_env %||% NULL
  if (is.environment(trace_env)) {
    trace_name <- scenario_spec$trace_name %||% "rule_trace"
    after_stats <- out_dt[, c(
      list(TimeInn = unique(get(config$time_in)), n = .N),
      stats::setNames(
        as.list(unlist(lapply(config$exposures, function(exp_name) mean(get(exp_name), na.rm = TRUE)))),
        paste0(config$exposures, "_after")
      )
    )]
    trace_dt <- merge(before_stats, after_stats, by = c("TimeInn", "n"), all = TRUE, sort = TRUE)
    for (exp_name in config$exposures) {
      before_col <- paste0(exp_name, "_before")
      after_col <- paste0(exp_name, "_after")
      ratio_col <- paste0(exp_name, "_ratio")
      trace_dt[, (ratio_col) := get(after_col) / get(before_col)]
    }
    if (exists(trace_name, envir = trace_env, inherits = FALSE)) {
      existing <- get(trace_name, envir = trace_env, inherits = FALSE)
      trace_dt <- data.table::rbindlist(list(existing, trace_dt), use.names = TRUE, fill = TRUE)
    }
    assign(trace_name, trace_dt[], envir = trace_env)
  }
  dt
}

prepare_counterfactual_dt <- function(mc_data, config) {
  dt <- data.table::as.data.table(data.table::copy(mc_data))
  force_levels_dt(dt, dt, intersect(config$factor_vars, names(dt)))
  data.table::setorderv(dt, c("uid", config$time_in))
  dt
}

simulate_longitudinal_panel <- function(dt,
                                        config,
                                        model_fit,
                                        exposure_mode = c("fixed", "modeled", "observed"),
                                        fixed_values = NULL,
                                        scenario_name = NULL,
                                        progress_fn = NULL,
                                        progress_label = NULL) {
  exposure_mode <- match.arg(exposure_mode)
  scenario_spec <- NULL
  if (identical(exposure_mode, "fixed")) {
    scenario_spec <- normalize_scenario(fixed_values, config, scenario_name = scenario_name)
    if (identical(scenario_spec$type, "static")) {
      vals <- scenario_spec$values
      for (i in seq_along(config$exposures)) {
        exp_name <- config$exposures[i]
        lag_name <- config$exposure_lags[i]
        dt[, (exp_name) := vals[[exp_name]]]
        if (lag_name %in% names(dt)) dt[, (lag_name) := vals[[exp_name]]]
      }
    }
  } else if (identical(exposure_mode, "modeled") && is.null(model_fit$exposure_models)) {
    stop("Modeled exposure simulation requires exposure_formulas in the config.")
  }

  base_time <- min(dt[[config$time_in]], na.rm = TRUE)
  times <- sort(unique(dt[[config$time_in]]))
  history_plan <- get_history_plan(config)

  prev_tvc_cols <- c("uid", config$tvc)
  prev_tvc <- dt[get(config$time_in) == base_time, ..prev_tvc_cols]
  data.table::setnames(prev_tvc, config$tvc, paste0(config$tvc, "_prev"))
  data.table::setkeyv(prev_tvc, "uid")

  prev_exp_cols <- c("uid", config$exposures)
  prev_exp <- dt[get(config$time_in) == base_time, ..prev_exp_cols]
  data.table::setnames(prev_exp, config$exposures, paste0(config$exposures, "_prev"))
  data.table::setkeyv(prev_exp, "uid")
  for (tt in times) {
    idx <- which(dt[[config$time_in]] == tt)
    id_example <- dt$uid[idx][1]
    if (tt > base_time) {
      dt <- populate_history_for_time(dt, config, history_plan, tt, times)
      dt <- carry_forward_passive_tvcs(dt, config, tt, prev_tvc = prev_tvc)

      for (v in config$tvc) {
        if (v %in% (config$passive_tvc_vars %||% character())) {
          next
        }
        dt[idx, (v) := simulate_from_model(
          model = model_fit$tvc_models[[v]],
          newdata = dt[idx],
          model_type = config$covtypes[[v]],
          support = model_fit$tvc_support[[v]],
          value_name = v,
          time_value = tt,
          id_value = id_example,
          categorical_sim = config$categorical_sim %||% "stochastic"
        )]
        assert_simulated_values(dt[[v]][idx], v, tt, id_value = id_example)
      }

      if (identical(exposure_mode, "modeled") ||
          (identical(exposure_mode, "fixed") && !is.null(scenario_spec) && identical(scenario_spec$type, "rule"))) {
        if (is.null(model_fit$exposure_models)) {
          stop("Rule-based exposure scenarios require fitted exposure models.")
        }
        for (v in config$exposures) {
          dt[idx, (v) := simulate_from_model(
            model = model_fit$exposure_models[[v]],
            newdata = dt[idx],
            model_type = config$exposure_types[[v]],
            support = model_fit$exposure_support[[v]],
            value_name = v,
            time_value = tt,
            id_value = id_example,
            categorical_sim = config$categorical_sim %||% "stochastic"
          )]
          assert_simulated_values(dt[[v]][idx], v, tt, id_value = id_example)
        }
      }
    }
    if (identical(exposure_mode, "fixed") && !is.null(scenario_spec) && identical(scenario_spec$type, "lookup")) {
      dt <- apply_lookup_scenario_at_time(dt, scenario_spec, config, tt)
    } else if (identical(exposure_mode, "fixed") && !is.null(scenario_spec) && identical(scenario_spec$type, "rule")) {
      dt <- apply_rule_scenario_at_time(dt, scenario_spec, config, tt)
    }

    dt[idx, hazard := safe_predict_response(
      model_fit$outcome_model,
      newdata = dt[idx],
      varname = "outcome hazard",
      time_value = tt,
      id_value = id_example
    )]
    dt[idx, hazard := pmin(pmax(hazard, 0), 1)]

    prev_tvc <- dt[get(config$time_in) == tt, ..prev_tvc_cols]
    data.table::setnames(prev_tvc, config$tvc, paste0(config$tvc, "_prev"))
    data.table::setkeyv(prev_tvc, "uid")

    prev_exp <- dt[get(config$time_in) == tt, ..prev_exp_cols]
    data.table::setnames(prev_exp, config$exposures, paste0(config$exposures, "_prev"))
    data.table::setkeyv(prev_exp, "uid")
    if (is.function(progress_fn)) progress_fn(sprintf("%s time %d/%d", progress_label, match(tt, times), length(times)))
  }

  dt[, cumrisk := 1 - cumprod(1 - hazard), by = uid]
  dt[, survival := 1 - cumrisk]
  if (!is.null(scenario_name)) dt[, scenario := scenario_name]
  extra_cols <- c("hazard", "cumrisk", "survival", "scenario", "joint_effect")
  trim_simulation_columns(dt, config, extra_cols = extra_cols)
}

simulate_modeled_exposome_natural_course <- function(mc_data, config, model_fit, progress_fn = NULL, progress_label = NULL) {
  dt <- prepare_counterfactual_dt(mc_data, config)
  simulate_longitudinal_panel(dt, config, model_fit, exposure_mode = "modeled", scenario_name = "natural_modeled_exposome", progress_fn = progress_fn, progress_label = progress_label)
}

simulate_observed_exposome_natural_course <- function(mc_data, config, model_fit, progress_fn = NULL, progress_label = NULL) {
  dt <- prepare_counterfactual_dt(mc_data, config)
  simulate_longitudinal_panel(dt, config, model_fit, exposure_mode = "observed", scenario_name = "natural_observed_exposome", progress_fn = progress_fn, progress_label = progress_label)
}

simulate_tvcqgcomp_scenario <- function(mc_data,
                                        config,
                                        model_fit,
                                        scenario_values,
                                        scenario_name = NULL,
                                        progress_fn = NULL,
                                        progress_label = NULL) {
  dt <- prepare_counterfactual_dt(mc_data, config)
  scenario_spec <- normalize_scenario(scenario_values, config, scenario_name = scenario_name)
  dt <- simulate_longitudinal_panel(dt, config, model_fit, exposure_mode = "fixed", fixed_values = scenario_values, scenario_name = scenario_name %||% "scenario", progress_fn = progress_fn, progress_label = progress_label)
  if (identical(scenario_spec$type, "static")) {
    vals <- scenario_spec$values
    dt[, joint_effect := vals[[config$exposures[1L]]]]
    for (exp_name in names(vals)) dt[, (paste0(exp_name, "_assigned")) := vals[[exp_name]]]
  } else {
    dt[, joint_effect := NA_real_]
    for (exp_name in config$exposures) dt[, (paste0(exp_name, "_assigned")) := get(exp_name)]
  }
  extra_cols <- c("hazard", "cumrisk", "survival", "scenario", "joint_effect")
  trim_simulation_columns(dt, config, extra_cols = extra_cols)
}

should_fit_meta_model <- function(exposure_scale, fit_meta = NULL) {
  if (!is.null(fit_meta)) {
    return(isTRUE(fit_meta))
  }
  identical(exposure_scale, "quantized")
}

validate_meta_intervention_support <- function(intervention_scenarios,
                                               config,
                                               fit_meta,
                                               call_name = "run_tvcqgcomp") {
  if (!isTRUE(fit_meta) || is.null(intervention_scenarios) || !length(intervention_scenarios)) {
    return(invisible(NULL))
  }

  unsupported <- names(intervention_scenarios)[vapply(
    names(intervention_scenarios),
    function(scn) {
      vals <- normalize_scenario(intervention_scenarios[[scn]], config, scenario_name = scn)
      identical(vals$type, "lookup") || identical(vals$type, "rule")
    },
    logical(1L)
  )]

  if (length(unsupported)) {
    stop(
      call_name,
      " does not support fit_meta = TRUE for lookup/rule interventions because no scalar joint_effect is defined for those scenarios. ",
      "Set fit_meta = FALSE and use risk_trajectory instead. Unsupported scenario(s): ",
      paste(unsupported, collapse = ", ")
    )
  }

  invisible(NULL)
}

should_fit_exposure_models <- function(config, natural_course, intervention_scenarios, fit_exposure_models = NULL) {
  if (!is.null(fit_exposure_models)) {
    return(isTRUE(fit_exposure_models))
  }
  natural_course <- normalize_natural_course_modes(natural_course)
  if (!identical(natural_course, "none")) {
    return(TRUE)
  }
  if (is.null(intervention_scenarios) || !length(intervention_scenarios)) {
    return(TRUE)
  }
  all(vapply(intervention_scenarios, function(vals) {
    if (is.list(vals) && identical(vals$type %||% NULL, "rule")) {
      return(FALSE)
    }
    all(config$exposures %in% names(vals))
  }, logical(1)))
}

make_final_risk_summary <- function(intervention_summary,
                                    natural_risk,
                                    natural_label = "natural") {
  intervention_summary <- data.table::as.data.table(data.table::copy(intervention_summary))

  has_natural <- !is.null(natural_risk) && is.finite(natural_risk) && !is.na(natural_risk)
  natural_summary <- NULL
  if (has_natural) {
    natural_summary <- data.table::data.table(
      scenario = natural_label,
      mean_final_risk = natural_risk,
      mean_final_risk_per_1000 = 1000 * natural_risk,
      mean_final_survival = 1 - natural_risk
    )
  }

  if (!is.null(natural_summary)) {
    data.table::rbindlist(list(natural_summary, intervention_summary), fill = TRUE)
  } else {
    intervention_summary
  }
}

flatten_final_risk_summary <- function(final_risk_summary) {
  dt <- data.table::as.data.table(data.table::copy(final_risk_summary))
  keep <- intersect(c("scenario", "mean_final_risk"), names(dt))
  dt <- dt[, ..keep]
  out <- list()
  for (i in seq_len(nrow(dt))) {
    scn <- dt$scenario[[i]]
    if ("mean_final_risk" %in% names(dt)) out[[paste0("risk_", scn)]] <- dt$mean_final_risk[[i]]
  }
  data.table::as.data.table(out)
}

summarize_risk_trajectory <- function(data, config, scenario_col = NULL, risk_col = "cumrisk") {
  dt <- data.table::as.data.table(data.table::copy(data))
  if (!(risk_col %in% names(dt))) stop("risk_col not found in data.")
  if (!"survival" %in% names(dt)) {
    dt[, survival := 1 - get(risk_col)]
  }
  grouping_columns <- c(
    if (!is.null(scenario_col) && scenario_col %in% names(dt)) scenario_col,
    config$time_in,
    if (!is.null(config$time_out) && config$time_out %in% names(dt)) config$time_out
  )
  dt[, .(
    mean_risk = mean(get(risk_col), na.rm = TRUE),
    mean_risk_per_1000 = 1000 * mean(get(risk_col), na.rm = TRUE),
    mean_survival = mean(survival, na.rm = TRUE),
    mean_survival_per_1000 = 1000 * mean(survival, na.rm = TRUE)
  ), by = grouping_columns]
}

flatten_risk_trajectory_endpoint <- function(risk_trajectory, config) {
  dt <- data.table::as.data.table(data.table::copy(risk_trajectory))
  time_col <- if (!is.null(config$time_out) && config$time_out %in% names(dt)) {
    config$time_out
  } else {
    config$time_in
  }
  if (is.null(time_col) || !(time_col %in% names(dt))) {
    stop("risk_trajectory does not contain the configured follow-up time column.")
  }
  if (!all(c("scenario", "mean_risk") %in% names(dt))) {
    stop("risk_trajectory must contain scenario and mean_risk columns.")
  }

  final_time <- max(dt[[time_col]], na.rm = TRUE)
  endpoint <- dt[get(time_col) == final_time, .(
    mean_risk = mean(mean_risk, na.rm = TRUE)
  ), by = scenario]

  values <- stats::setNames(
    as.list(endpoint$mean_risk),
    paste0("risk_", endpoint$scenario)
  )
  data.table::as.data.table(values)
}

as_risk_curve_df <- function(run_obj,
                             include_natural = TRUE,
                             natural_label = "Natural Course",
                             time_col = NULL,
                             scenario_col = "scenario") {
  if (!is.list(run_obj)) {
    stop("run_obj must be a tvcQGcomp run object.")
  }
  risk_trajectory <- run_obj$risk_trajectory %||% run_obj$intervention_curve
  if (is.null(risk_trajectory) && is.null(run_obj$natural_curve)) {
    stop("run_obj does not contain curve summaries.")
  }

  config <- run_obj$config %||% NULL
  if (is.null(time_col)) {
    time_col <- config$time_in %||% "TimeInn"
  }

  out <- list()
  if (isTRUE(include_natural) && !is.null(run_obj$natural_curve)) {
    nat <- data.table::as.data.table(data.table::copy(run_obj$natural_curve))
    nat[, (scenario_col) := natural_label]
    out[[length(out) + 1L]] <- nat
  }
  if (!is.null(risk_trajectory)) {
    intv <- data.table::as.data.table(data.table::copy(risk_trajectory))
    if (!(scenario_col %in% names(intv))) {
      stop("scenario_col not found in risk_trajectory.")
    }
    out[[length(out) + 1L]] <- intv
  }
  if (!length(out)) {
    stop("No risk-curve data available.")
  }

  curve_dt <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  if (!(time_col %in% names(curve_dt))) {
    stop("time_col not found in curve data.")
  }
  data.table::setorderv(curve_dt, c(scenario_col, time_col))
  curve_dt[]
}


make_boot_curve_result <- function(run_obj,
                                   boot,
                                   include_natural = TRUE,
                                   natural_label = "Natural Course",
                                   scenario_col = "scenario",
                                   time_col = NULL) {
  curve_dt <- as_risk_curve_df(
    run_obj,
    include_natural = include_natural,
    natural_label = natural_label,
    time_col = time_col,
    scenario_col = scenario_col
  )
  curve_dt[, boot := boot]
  data.table::setcolorder(curve_dt, c("boot", setdiff(names(curve_dt), "boot")))
  curve_dt[]
}

summarize_boot_risk_curves <- function(boot_res,
                                       scenario_col = "scenario",
                                       time_col = NULL,
                                       conf_level = 0.95) {
  br <- if (is.list(boot_res) && !is.null(boot_res$curve_results)) {
    boot_res$curve_results
  } else {
    NULL
  }
  if (is.null(br)) {
    stop("boot_res does not contain stored curve_results. Re-run the bootstrap with keep_boot_curves = TRUE.")
  }
  dt <- data.table::as.data.table(data.table::copy(br))
  if (is.null(time_col)) {
    time_col <- if ("TimeInn" %in% names(dt)) "TimeInn" else if ("time" %in% names(dt)) "time" else NULL
  }
  if (is.null(time_col) || !(time_col %in% names(dt))) {
    stop("time_col not found in curve_results.")
  }
  if (!(scenario_col %in% names(dt))) {
    stop("scenario_col not found in curve_results.")
  }

  alpha <- 1 - conf_level
  lower_p <- alpha / 2
  upper_p <- 1 - alpha / 2
  metric_cols <- intersect(
    c("mean_risk", "mean_survival"),
    names(dt)
  )
  if (!length(metric_cols)) {
    stop("No curve metrics found in curve_results.")
  }

  stats_list <- lapply(metric_cols, function(metric) {
    dt[, .(
      metric = metric,
      estimate = mean(get(metric), na.rm = TRUE),
      median = stats::median(get(metric), na.rm = TRUE),
      std_error = stats::sd(get(metric), na.rm = TRUE),
      conf.low = stats::quantile(get(metric), lower_p, na.rm = TRUE),
      conf.high = stats::quantile(get(metric), upper_p, na.rm = TRUE)
    ), by = c(scenario_col, time_col)]
  })
  out <- data.table::rbindlist(stats_list, use.names = TRUE, fill = TRUE)
  data.table::setorderv(out, c(scenario_col, time_col, "metric"))
  out[]
}

plot_cumulative_risk_trajectory <- function(x,
                                            include_natural = TRUE,
                                            natural_label = "Natural Course",
                                            include_average = FALSE,
                                            average_scenarios = paste0("Q", 0:3),
                                            average_label = "Average (Q0-Q3)",
                                            average_col = "black",
                                            scenario_col = "scenario",
                                            time_col = "TimeInn",
                                            risk_col = "mean_risk",
                                            percent = TRUE,
                                            scenario_order = NULL,
                                            scenario_labels = NULL,
                                            cols = NULL,
                                            include_points = TRUE,
                                            pch = 16,
                                            lwd = 2,
                                            xlab = "Follow-up time",
                                            ylab = if (percent) "Cumulative Risk (pct)" else "Cumulative Risk",
                                            main = "Cumulative risk over time",
                                            legend_pos = "topleft",
                                            ...) {
  dt <- NULL
  if (is.list(x) && !data.table::is.data.table(x)) {
    dt <- as_risk_curve_df(
      x,
      include_natural = include_natural,
      natural_label = natural_label,
      time_col = time_col,
      scenario_col = scenario_col
    )
  } else {
    dt <- data.table::as.data.table(data.table::copy(x))
  }

  if (!(time_col %in% names(dt))) stop("time_col not found in data.")
  if (!(risk_col %in% names(dt))) stop("risk_col not found in data.")
  if (!(scenario_col %in% names(dt))) dt[, (scenario_col) := "risk"]

  dt <- dt[!is.na(get(time_col)) & !is.na(get(risk_col))]
  if (!nrow(dt)) stop("No non-missing rows available to plot.")

  if (isTRUE(percent)) {
    dt[, plot_risk := 100 * get(risk_col)]
  } else {
    dt[, plot_risk := get(risk_col)]
  }

  if (isTRUE(include_average)) {
    average_scenarios <- as.character(average_scenarios)
    present_average <- intersect(average_scenarios, unique(as.character(dt[[scenario_col]])))
    if (!length(present_average)) {
      warning("No average_scenarios were found in the trajectory data; average line was not added.", call. = FALSE)
    } else {
      if (length(present_average) < length(average_scenarios)) {
        warning(
          "Only these average_scenarios were found: ",
          paste(present_average, collapse = ", "),
          call. = FALSE
        )
      }
      avg_dt <- dt[as.character(get(scenario_col)) %in% present_average,
        .(plot_risk = mean(plot_risk, na.rm = TRUE)),
        by = time_col
      ]
      avg_dt[, (scenario_col) := average_label]
      avg_dt[, (risk_col) := if (isTRUE(percent)) plot_risk / 100 else plot_risk]
      dt <- data.table::rbindlist(list(dt, avg_dt), use.names = TRUE, fill = TRUE)
    }
  }

  scenarios <- unique(as.character(dt[[scenario_col]]))
  if (!is.null(scenario_order)) {
    scenario_order <- as.character(scenario_order)
    scenarios <- c(intersect(scenario_order, scenarios), setdiff(scenarios, scenario_order))
  }

  if (is.null(cols)) {
    palette_cols <- grDevices::hcl.colors(max(length(scenarios), 3L), palette = "Dark 3")
    cols <- stats::setNames(palette_cols[seq_along(scenarios)], scenarios)
    if (natural_label %in% scenarios) cols[[natural_label]] <- "#E41A1C"
    if (average_label %in% scenarios) cols[[average_label]] <- average_col
  } else if (is.null(names(cols))) {
    if (length(cols) < length(scenarios)) stop("Unnamed cols must have length at least equal to the number of scenarios.")
    cols <- stats::setNames(cols[seq_along(scenarios)], scenarios)
  } else {
    miss_cols <- setdiff(scenarios, names(cols))
    if (length(miss_cols)) {
      add_cols <- grDevices::hcl.colors(length(miss_cols), palette = "Dark 3")
      cols <- c(cols, stats::setNames(add_cols, miss_cols))
    }
    cols <- cols[scenarios]
  }

  legend_labels <- scenarios
  if (!is.null(scenario_labels)) {
    if (is.null(names(scenario_labels))) {
      if (length(scenario_labels) != length(scenarios)) stop("Unnamed scenario_labels must have length equal to the number of scenarios.")
      legend_labels <- as.character(scenario_labels)
    } else {
      idx <- match(scenarios, names(scenario_labels))
      replace <- !is.na(idx)
      legend_labels[replace] <- as.character(scenario_labels[idx[replace]])
    }
  }
  ylim <- range(dt$plot_risk, na.rm = TRUE)
  if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)

  graphics::plot(NA,
    xlim = range(dt[[time_col]], na.rm = TRUE),
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    yaxt = if (percent) "n" else "s",
    ...
  )
  if (isTRUE(percent)) {
    yticks <- graphics::axTicks(2)
    graphics::axis(2, at = yticks, labels = paste0(format(round(yticks, 1), nsmall = 1), "%"))
  }

  for (scn in scenarios) {
    tmp <- dt[get(scenario_col) == scn]
    data.table::setorderv(tmp, time_col)
    graphics::lines(tmp[[time_col]], tmp$plot_risk, col = cols[[scn]], lwd = lwd)
    if (isTRUE(include_points)) {
      graphics::points(tmp[[time_col]], tmp$plot_risk, col = cols[[scn]], pch = pch)
    }
  }

  graphics::legend(
    legend_pos,
    legend = legend_labels,
    col = unname(cols[scenarios]),
    lwd = lwd,
    pch = if (isTRUE(include_points)) pch else NA,
    bty = "n"
  )
  invisible(dt[])
}

plot_survival_trajectory <- function(x,
                                     include_natural = TRUE,
                                     natural_label = "Natural Course",
                                     include_average = FALSE,
                                     average_scenarios = paste0("Q", 0:3),
                                     average_label = "Average (Q0-Q3)",
                                     average_col = "black",
                                     scenario_col = "scenario",
                                     time_col = "TimeInn",
                                     survival_col = "mean_survival",
                                     percent = TRUE,
                                     scenario_order = NULL,
                                     scenario_labels = NULL,
                                     cols = NULL,
                                     include_points = TRUE,
                                     pch = 16,
                                     lwd = 2,
                                     xlab = "Follow-up time",
                                     ylab = if (percent) "Survival (%)" else "Survival probability",
                                     main = "Survival over time",
                                     legend_pos = "bottomleft",
                                     ...) {
  plot_cumulative_risk_trajectory(
    x = x,
    include_natural = include_natural,
    natural_label = natural_label,
    include_average = include_average,
    average_scenarios = average_scenarios,
    average_label = average_label,
    average_col = average_col,
    scenario_col = scenario_col,
    time_col = time_col,
    risk_col = survival_col,
    percent = percent,
    scenario_order = scenario_order,
    scenario_labels = scenario_labels,
    cols = cols,
    include_points = include_points,
    pch = pch,
    lwd = lwd,
    xlab = xlab,
    ylab = ylab,
    main = main,
    legend_pos = legend_pos,
    ...
  )
}

#' @export
plot.tvcqgcomp_run <- function(x, ...) {
  plot_cumulative_risk_trajectory(x, ...)
}

normalize_meta_targets <- function(x) {
  unique(unlist(lapply(x, function(target) {
    target <- match.arg(target, c("HR", "HD", "RR", "RD", "hazard", "cumrisk_final", "both"))
    if (identical(target, "hazard")) return("HR")
    if (identical(target, "cumrisk_final")) return("RR")
    if (identical(target, "both")) return(c("HR", "RR"))
    target
  }), use.names = FALSE))
}

resolve_meta_specs <- function(spec, targets, arg_name = "meta specification") {
  if (length(targets) == 1L) {
    out <- list(spec)
    names(out) <- targets
    return(out)
  }
  if (is.null(spec)) {
    out <- vector("list", length(targets))
    names(out) <- targets
    return(out)
  }
  if (!is.list(spec)) {
    stop("When multiple meta targets are requested, ", arg_name, " must be NULL or a list named by target.")
  }
  if (is.null(names(spec))) {
    if (length(spec) != length(targets)) {
      stop("Unnamed ", arg_name, " list must match the number of meta targets.")
    }
    names(spec) <- targets
    return(spec)
  }
  if (!all(targets %in% names(spec))) {
    stop("Missing ", arg_name, " entries for: ", paste(setdiff(targets, names(spec)), collapse = ", ") )
  }
  spec[targets]
}

normalize_natural_course_modes <- function(x) {
  unique(unlist(lapply(x, function(mode) {
    mode <- match.arg(mode, c("none", "calibration", "observed", "observed_exposome", "modeled_exposome", "modeled"))
    if (identical(mode, "observed")) return("calibration")
    if (identical(mode, "modeled")) return("modeled_exposome")
    mode
  }), use.names = FALSE))
}

make_default_meta_formula <- function(config, meta_target) {
  if (identical(meta_target, "HR") || identical(meta_target, "HD")) {
    return(stats::as.formula(paste0("hazard ~ joint_effect + factor(", config$time_in, ")")))
  }
  stats::as.formula("cumrisk ~ joint_effect")
}

fit_pooled_meta_model <- function(simulated_data,
                                  config,
                                  meta_formula = NULL,
                                  meta_family = NULL,
                                  meta_target = config$meta_target %||% "HR") {
  dt <- data.table::as.data.table(data.table::copy(simulated_data))
  meta_target <- normalize_meta_targets(meta_target)
  if (length(meta_target) != 1L) {
    stop("fit_pooled_meta_model expects a single meta target.")
  }
  meta_target <- meta_target[[1L]]

  if (identical(meta_target, "HR") || identical(meta_target, "HD")) {
    if (is.null(meta_formula)) {
      meta_formula <- make_default_meta_formula(config, meta_target)
    }
    if (is.null(meta_family)) {
      meta_family <- if (identical(meta_target, "HD")) {
        stats::gaussian("identity")
      } else {
        stats::quasibinomial("logit")
      }
    }
    if (any(!is.finite(dt$hazard))) {
      stop("Simulated data contains invalid hazard values. Inspect the scenario outputs before fitting the meta model.")
    }
    fit <- stats::glm(meta_formula, data = dt, family = meta_family)
    attr(fit, "meta_target") <- meta_target
    return(fit)
  }

  final_time <- max(dt[[config$time_in]], na.rm = TRUE)
  final_dt <- dt[get(config$time_in) == final_time]
  if (is.null(meta_formula)) {
    meta_formula <- make_default_meta_formula(config, meta_target)
  }
  if (is.null(meta_family)) {
    meta_family <- if (identical(meta_target, "RD")) {
      stats::quasibinomial("identity")
    } else {
      stats::quasibinomial("log")
    }
  }
  if (any(!is.finite(final_dt$cumrisk))) {
    stop("Simulated data contains invalid cumrisk values at the final time. Inspect the scenario outputs before fitting the meta model.")
  }
  fit <- stats::glm(meta_formula, data = final_dt, family = meta_family)
  attr(fit, "meta_target") <- meta_target
  fit
}

fit_meta_models <- function(simulated_data,
                            config,
                            meta_formula = NULL,
                            meta_family = NULL,
                            meta_target = config$meta_target %||% "HR") {
  targets <- normalize_meta_targets(meta_target)
  formula_specs <- resolve_meta_specs(meta_formula, targets, arg_name = "meta_formula")
  family_specs <- resolve_meta_specs(meta_family, targets, arg_name = "meta_family")
  models <- setNames(vector("list", length(targets)), targets)
  for (target in targets) {
    models[[target]] <- fit_pooled_meta_model(
      simulated_data = simulated_data,
      config = config,
      meta_formula = formula_specs[[target]],
      meta_family = family_specs[[target]],
      meta_target = target
    )
  }
  models
}

summarize_meta_model_effects <- function(meta_models) {
  if (is.null(meta_models) || !length(meta_models)) {
    return(NULL)
  }

  ratio_targets <- c("HR", "RR")
  rows <- lapply(names(meta_models), function(target) {
    model <- meta_models[[target]]
    coefs <- stats::coef(model)
    link_estimate <- if ("joint_effect" %in% names(coefs)) {
      unname(coefs[["joint_effect"]])
    } else {
      NA_real_
    }
    effect_estimate <- if (target %in% ratio_targets) {
      exp(link_estimate)
    } else {
      link_estimate
    }

    data.table::data.table(
      meta_target = target,
      estimate = effect_estimate,
      increment = 1
    )
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

run_tvcqgcomp <- function(data,
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
                          workshop_contrasts = FALSE,
                          workshop_reference = NULL,
                          workshop_scenarios = NULL,
                          workshop_include_natural = FALSE,
                          verbose = FALSE) {
  meta_formula <- if (is.list(meta_formula)) with_formula_helpers_list(meta_formula) else with_formula_helpers(meta_formula)
  natural_course <- normalize_natural_course_modes(natural_course)
  if (is.null(intervention_scenarios)) intervention_scenarios <- make_q_scenarios(config)
  fit_meta <- should_fit_meta_model(exposure_scale, fit_meta)
  validate_meta_intervention_support(intervention_scenarios, config, fit_meta, call_name = "run_tvcqgcomp")
  fit_exposure_models <- should_fit_exposure_models(config, natural_course, intervention_scenarios, fit_exposure_models)

  format_progress_time <- function(seconds) {
    if (is.character(seconds)) return(seconds)
    if (!is.finite(seconds) || is.na(seconds)) return("estimating")
    seconds <- max(0, round(seconds))
    hrs <- seconds %/% 3600
    mins <- (seconds %% 3600) %/% 60
    secs <- seconds %% 60
    if (hrs > 0) sprintf("%02d:%02d:%02d", hrs, mins, secs) else sprintf("%02d:%02d", mins, secs)
  }

  print_progress_bar <- function(done, total, start_time, label = "tvcQGcomp") {
    total <- max(1L, as.integer(total))
    done <- min(max(0L, as.integer(done)), total)
    width <- 24L
    frac <- done / total
    filled <- as.integer(floor(width * frac))
    empty <- width - filled
    elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    eta_sec <- if (done > 0L) elapsed_sec / done * (total - done) else "estimating"
    green_on <- "\033[32m"
    color_off <- "\033[39m"
    bar <- paste0(green_on, strrep("=", filled), color_off, strrep("-", empty))
    pct <- sprintf("%3d%%", round(frac * 100))
    cat(sprintf("%s [%s] %s %d/%d | elapsed %s | ETA %s\n",
                label,
                bar,
                pct,
                done,
                total,
                format_progress_time(elapsed_sec),
                format_progress_time(eta_sec)))
  }

  prep <- prepare_exposure_scale_data(data, config, exposure_scale = exposure_scale, q = config$q_levels, breaks = quantization_breaks)
  dt <- prepare_tvcqgcomp_data(prep$data, config)
  force_levels_dt(dt, dt, intersect(config$factor_vars, names(dt)))

  n_times <- length(unique(dt[[config$time_in]]))
  n_passive_tvc <- length(config$passive_tvc_vars %||% character())
  n_tvc_models <- length(config$tvc) - n_passive_tvc
  n_tvc_steps <- n_tvc_models + n_passive_tvc
  n_exposure_models <- if (isTRUE(fit_exposure_models) && !is.null(config$exposure_formulas)) length(config$exposures) else 0L
  n_meta_models <- if (isTRUE(fit_meta)) length(normalize_meta_targets(meta_target)) else 0L
  progress_total <- 1L + 1L + n_tvc_steps + n_exposure_models + 1L + n_meta_models + 1L
  if ("calibration" %in% natural_course) progress_total <- progress_total + 1L
  if ("observed_exposome" %in% natural_course) progress_total <- progress_total + n_times
  if ("modeled_exposome" %in% natural_course) progress_total <- progress_total + n_times
  progress_total <- progress_total + length(names(intervention_scenarios)) * n_times

  start_time <- Sys.time()
  progress_done <- 0L
  progress_tick <- function(label = "tvcQGcomp") {
    progress_done <<- progress_done + 1L
    if (isTRUE(verbose)) print_progress_bar(progress_done, progress_total, start_time, label = label)
  }

  if (isTRUE(verbose)) {
    plan_parts <- c(
      "prepare data",
      "fit 1 outcome model",
      sprintf("fit %d confounder model(s)", n_tvc_models)
    )
    if (n_passive_tvc > 0L) {
      plan_parts <- c(plan_parts, sprintf("retain %d passive confounder path(s)", n_passive_tvc))
    }
    plan_parts <- c(
      plan_parts,
      sprintf("fit %d exposure model(s)", n_exposure_models),
      "build MC skeleton",
      sprintf(
        "simulate %d natural step(s)",
        (if ("calibration" %in% natural_course) 1L else 0L) +
          (if ("observed_exposome" %in% natural_course) n_times else 0L) +
          (if ("modeled_exposome" %in% natural_course) n_times else 0L)
      ),
      sprintf("simulate %d scenario step(s)", length(names(intervention_scenarios)) * n_times),
      sprintf("fit %d meta model(s)", n_meta_models),
      "finalize result"
    )
    cat("tvcQGcomp plan: ", paste(plan_parts, collapse = " + "), "\n", sep = "")
    print_progress_bar(0L, progress_total, start_time, label = "starting")
  }

  progress_tick("prepare data")
  model_fit <- fit_tvcqgcomp_models(dt, config, progress_fn = progress_tick, fit_exposure_models = fit_exposure_models)
  mc_data <- create_mc_skeleton(dt, config, mc_size = mc_size, seed = seed, replace = replace_mc)
  progress_tick("build MC skeleton")

  natural_calibration <- NULL
  natural_observed_exposome <- NULL
  natural_modeled_exposome <- NULL
  if ("calibration" %in% natural_course) {
    natural_calibration <- predict_natural_course(dt, config, model_fit, mode = "calibration")
    progress_tick("simulate natural calibration")
  }
  if ("observed_exposome" %in% natural_course) {
    natural_observed_exposome <- simulate_observed_exposome_natural_course(mc_data, config, model_fit, progress_fn = progress_tick, progress_label = "simulate natural observed exposome")
  }
  if ("modeled_exposome" %in% natural_course) {
    natural_modeled_exposome <- simulate_modeled_exposome_natural_course(mc_data, config, model_fit, progress_fn = progress_tick, progress_label = "simulate natural modeled exposome")
  }

  natural_default <- if ("observed_exposome" %in% natural_course) {
    natural_observed_exposome
  } else if ("calibration" %in% natural_course) {
    natural_calibration
  } else {
    natural_modeled_exposome
  }

  scenario_list <- lapply(names(intervention_scenarios), function(scn) {
    vals <- normalize_scenario(intervention_scenarios[[scn]], config, scenario_name = scn)
    sim_dt <- simulate_tvcqgcomp_scenario(mc_data, config, model_fit, vals, scenario_name = scn, progress_fn = progress_tick, progress_label = paste0("simulate scenario ", scn))
    joint_value <- if (identical(vals$type, "lookup") || identical(vals$type, "rule")) {
      NA_real_
    } else if (is.null(joint_effect_fn)) {
      vals$values[[config$exposures[1L]]]
    } else {
      joint_effect_fn(vals$values)
    }
    sim_dt[, joint_effect := joint_value]
    sim_dt
  })
  names(scenario_list) <- names(intervention_scenarios)

  intervention_data <- data.table::rbindlist(scenario_list, fill = TRUE)

  workshop_summary <- NULL

  meta_models <- NULL
  meta_targets <- character()
  meta_model <- NULL

  if (isTRUE(fit_meta)) {
    meta_models <- fit_meta_models(
      intervention_data,
      config,
      meta_formula = meta_formula,
      meta_family = meta_family,
      meta_target = meta_target
    )
    for (target in names(meta_models)) {
      progress_tick(paste0("fit meta model: ", target))
    }
    meta_targets <- names(meta_models %||% list())
    meta_model <- if (length(meta_targets)) meta_models[[meta_targets[[1L]]]] else NULL
  }

  meta_effect_summary <- summarize_meta_model_effects(meta_models)

  out <- list(
    natural = natural_default,
      calibration = list(
        calibration = natural_calibration,
        calibration_curve = if (is.null(natural_calibration)) NULL else summarize_risk_trajectory(natural_calibration, config),
        observed = natural_calibration,
        observed_curve = if (is.null(natural_calibration)) NULL else summarize_risk_trajectory(natural_calibration, config),
        observed_exposome = natural_observed_exposome,
        observed_exposome_curve = if (is.null(natural_observed_exposome)) NULL else summarize_risk_trajectory(natural_observed_exposome, config),
        modeled_exposome = natural_modeled_exposome,
        modeled_exposome_curve = if (is.null(natural_modeled_exposome)) NULL else summarize_risk_trajectory(natural_modeled_exposome, config),
        modeled = natural_modeled_exposome,
        modeled_curve = if (is.null(natural_modeled_exposome)) NULL else summarize_risk_trajectory(natural_modeled_exposome, config)
      ),
    mc_data = mc_data,
    interventions = scenario_list,
    intervention_data = intervention_data,
    model_fit = model_fit,
    meta_model = meta_model,
    meta_models = meta_models,
    meta_effect_summary = meta_effect_summary,
    meta_target = meta_targets,
    natural_curve = if (is.null(natural_default)) NULL else summarize_risk_trajectory(natural_default, config),
    risk_trajectory = summarize_risk_trajectory(intervention_data, config, scenario_col = "scenario"),
    workshop_summary = workshop_summary,
    config = config,
    fit_meta = fit_meta,
    fit_exposure_models = fit_exposure_models
  )
  for (nm in names(scenario_list)) out[[nm]] <- scenario_list[[nm]]
  class(out) <- c("tvcqgcomp_run", "list")
  progress_tick("finalize result")
  out
}

validate_bootstrap_parallel_args <- function(parallel,
                                             n_workers,
                                             batch_size = NULL,
                                             audit_mode = FALSE,
                                             keep_boot_objects = FALSE,
                                             keep_boot_curves = FALSE,
                                             call_name = "bootstrap_tvcqgcomp") {
  if (!isTRUE(parallel)) {
    return(list(
      n_workers = as.integer(n_workers),
      batch_size = if (is.null(batch_size)) NULL else as.integer(batch_size)
    ))
  }
  if (length(n_workers) != 1L || is.na(n_workers) || !is.finite(n_workers) || n_workers < 1) {
    stop(call_name, " requires n_workers to be a single positive number when parallel = TRUE.")
  }
  n_workers <- as.integer(n_workers)
  if (n_workers < 1L) {
    stop(call_name, " requires n_workers >= 1 when parallel = TRUE.")
  }
  if (!is.null(batch_size)) {
    if (length(batch_size) != 1L || is.na(batch_size) || !is.finite(batch_size) || batch_size < 1) {
      stop(call_name, " requires batch_size to be NULL or a single positive number.")
    }
    batch_size <- as.integer(batch_size)
    if (batch_size < 1L) {
      stop(call_name, " requires batch_size >= 1 when supplied.")
    }
  }
  if (isTRUE(audit_mode) || isTRUE(keep_boot_objects) || isTRUE(keep_boot_curves)) {
    warning(
      call_name,
      " is running with parallel = TRUE and bootstrap object retention enabled. ",
      "This can use substantial RAM because each worker may return large per-bootstrap objects."
    )
  }
  list(n_workers = n_workers, batch_size = batch_size)
}

split_parallel_bootstrap_ids <- function(remaining_ids,
                                         n_workers,
                                         batch_size = NULL) {
  remaining_ids <- as.integer(remaining_ids)
  if (!length(remaining_ids)) {
    return(list())
  }
  if (is.null(batch_size)) {
    chunk_size <- length(remaining_ids)
  } else {
    chunk_size <- max(as.integer(n_workers), as.integer(batch_size))
  }
  split(remaining_ids, ceiling(seq_along(remaining_ids) / chunk_size))
}

run_single_bootstrap_worker <- function(b,
                                        dt,
                                        config,
                                        ids,
                                        level_vars,
                                        mc_size,
                                        intervention_scenarios,
                                        natural_course,
                                        seed,
                                        replace_mc,
                                        meta_formula,
                                        meta_family,
                                        meta_target,
                                        fit_meta,
                                        fit_exposure_models,
                                        exposure_scale,
                                        quantization_breaks,
                                        joint_effect_fn,
                                        keep_boot_objects,
                                        keep_boot_curves,
                                        dt_threads = NULL) {
  if (!is.null(dt_threads)) {
    old_dt_threads <- tryCatch(data.table::getDTthreads(), error = function(e) NA_integer_)
    data.table::setDTthreads(as.integer(dt_threads))
    if (!is.na(old_dt_threads)) {
      on.exit(data.table::setDTthreads(old_dt_threads), add = TRUE)
    }
  }
  set.seed(seed + b)
  boot_ids <- sample(ids, size = length(ids), replace = TRUE)
  boot_map <- data.table::data.table(tmp_id = boot_ids, boot_id = seq_along(boot_ids))
  data.table::setnames(boot_map, 'tmp_id', config$id)
  boot_dt <- merge(boot_map, dt, by = config$id, allow.cartesian = TRUE)
  boot_dt[, original_id := get(config$id)]
  boot_dt[, (config$id) := paste0(boot_id, '_', original_id)]
  boot_dt <- force_levels_dt(boot_dt, dt, level_vars)

  out <- tryCatch(
    run_tvcqgcomp(
      data = boot_dt,
      config = config,
      mc_size = mc_size,
      intervention_scenarios = intervention_scenarios,
      natural_course = natural_course,
      seed = seed + 100000L + b,
      replace_mc = replace_mc,
      meta_formula = meta_formula,
      meta_family = meta_family,
      meta_target = meta_target,
      fit_meta = fit_meta,
      fit_exposure_models = fit_exposure_models,
      exposure_scale = exposure_scale,
      quantization_breaks = quantization_breaks,
      joint_effect_fn = joint_effect_fn
    ),
    error = function(e) e
  )

  if (inherits(out, 'error')) {
    return(list(
      result = data.table::data.table(boot = b, status = paste('error:', conditionMessage(out))),
      object = if (keep_boot_objects) list(boot_data = boot_dt, error = conditionMessage(out)) else NULL,
      curves = NULL
    ))
  }

  if (isTRUE(fit_meta)) {
    meta_models <- out$meta_models %||% list(default = out$meta_model)
    if (length(meta_models) == 1L) {
      coefs <- stats::coef(meta_models[[1L]])
      coef_dt <- data.table::as.data.table(as.list(coefs))
      data.table::setnames(coef_dt, names(coef_dt), paste0('coef_', names(coef_dt)))
    } else {
      coef_vals <- unlist(lapply(names(meta_models), function(target) {
        coefs <- stats::coef(meta_models[[target]])
        stats::setNames(as.list(coefs), paste0('coef_', target, '_', names(coefs)))
      }), recursive = FALSE, use.names = TRUE)
      coef_dt <- data.table::as.data.table(coef_vals)
    }
  } else {
    coef_dt <- flatten_risk_trajectory_endpoint(out$risk_trajectory, config)
  }
  coef_dt[, boot := b]
  coef_dt[, status := 'success']
  data.table::setcolorder(coef_dt, c('boot', 'status', setdiff(names(coef_dt), c('boot', 'status'))))

  audit_object <- NULL
  if (keep_boot_objects) {
    audit_object <- list(
      boot_data = boot_dt,
      mc_data = out$mc_data,
      natural = out$natural,
      natural_calibration = out$calibration$calibration,
      natural_observed = out$calibration$observed,
      natural_modeled_exposome = out$calibration$modeled_exposome,
      natural_modeled = out$calibration$modeled,
      interventions = out$interventions,
      intervention_data = out$intervention_data,
      natural_curve = out$natural_curve,
      risk_trajectory = out$risk_trajectory,
      meta_model = out$meta_model,
      meta_models = out$meta_models
    )
  }

  curve_result <- NULL
  if (keep_boot_curves) {
    curve_result <- tryCatch(
      make_boot_curve_result(out, boot = b),
      error = function(e) NULL
    )
  }

  rm(out, boot_dt, boot_ids, boot_map)
  gc(verbose = FALSE)
  list(result = coef_dt[], object = audit_object, curves = curve_result)
}

bootstrap_tvcqgcomp <- function(data,
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
  dt <- prepare_tvcqgcomp_data(data, config)
  fit_meta <- should_fit_meta_model(exposure_scale, fit_meta)
  fit_exposure_models <- should_fit_exposure_models(config, natural_course, intervention_scenarios, fit_exposure_models)
  validate_meta_intervention_support(intervention_scenarios, config, fit_meta, call_name = "bootstrap_tvcqgcomp")
  ids <- unique(dt[[config$id]])
  level_vars <- intersect(config$factor_vars, names(dt))
  start_time <- Sys.time()
  keep_boot_objects <- isTRUE(keep_boot_objects) || isTRUE(audit_mode)
  keep_boot_curves <- isTRUE(keep_boot_curves) || keep_boot_objects
  parallel_args <- validate_bootstrap_parallel_args(
    parallel = parallel,
    n_workers = n_workers,
    batch_size = batch_size,
    audit_mode = audit_mode,
    keep_boot_objects = keep_boot_objects,
    keep_boot_curves = keep_boot_curves,
    call_name = "bootstrap_tvcqgcomp"
  )
  n_workers <- parallel_args$n_workers
  batch_size <- parallel_args$batch_size

  run_single_bootstrap <- run_single_bootstrap_worker

  results <- vector('list', n_boot)
  audit_objects <- if (keep_boot_objects) vector('list', n_boot) else NULL
  curve_results <- if (keep_boot_curves) vector('list', n_boot) else NULL

  restore_checkpoint <- function() {
    if (is.null(checkpoint_file) || !isTRUE(resume_from_checkpoint) || !file.exists(checkpoint_file)) {
      return(0L)
    }
    payload <- readRDS(checkpoint_file)
    if (is.list(payload) && !is.null(payload$results) && !data.table::is.data.table(payload$results)) {
      completed_n <- min(length(payload$results), n_boot)
      if (completed_n > 0L) {
        results[seq_len(completed_n)] <<- payload$results[seq_len(completed_n)]
      }
      if (keep_boot_objects && !is.null(payload$audit_objects)) {
        audit_n <- min(length(payload$audit_objects), n_boot)
        if (audit_n > 0L) {
          audit_objects[seq_len(audit_n)] <<- payload$audit_objects[seq_len(audit_n)]
        }
      }
      if (keep_boot_curves && !is.null(payload$curve_results)) {
        curve_n <- min(length(payload$curve_results), n_boot)
        if (curve_n > 0L) {
          curve_results[seq_len(curve_n)] <<- payload$curve_results[seq_len(curve_n)]
        }
      }
      return(min(as.integer(payload$last_completed %||% completed_n), n_boot))
    }
    if (is.list(payload) && data.table::is.data.table(payload$results)) {
      payload <- payload$results
    }
    if (data.table::is.data.table(payload) || is.data.frame(payload)) {
      dt_payload <- data.table::as.data.table(payload)
      if (!("boot" %in% names(dt_payload))) return(0L)
      boot_ids <- sort(unique(stats::na.omit(as.integer(dt_payload$boot))))
      for (b in boot_ids[boot_ids >= 1L & boot_ids <= n_boot]) {
        results[[b]] <<- dt_payload[boot == b]
      }
      return(if (length(boot_ids)) max(boot_ids) else 0L)
    }
    0L
  }

  save_checkpoint <- function(last_completed) {
    if (!is.null(checkpoint_file)) {
      completed <- seq_len(last_completed)
      payload <- list(results = results[completed], last_completed = last_completed)
      if (keep_boot_objects) payload$audit_objects <- audit_objects[completed]
      if (keep_boot_curves) payload$curve_results <- curve_results[completed]
      saveRDS(payload, checkpoint_file)
    }
  }
  checkpoint_completed <- restore_checkpoint()
  if (isTRUE(verbose) && checkpoint_completed > 0L) {
    cat('Resuming from checkpoint after bootstrap', checkpoint_completed, '\n')
  }
  next_boot <- checkpoint_completed + 1L
  format_progress_time <- function(seconds) {
    if (!is.finite(seconds) || is.na(seconds)) return('NA')
    seconds <- max(0, round(seconds))
    hrs <- seconds %/% 3600
    mins <- (seconds %% 3600) %/% 60
    secs <- seconds %% 60
    if (hrs > 0) {
      sprintf('%02d:%02d:%02d', hrs, mins, secs)
    } else {
      sprintf('%02d:%02d', mins, secs)
    }
  }

  print_progress_bar <- function(done, total, start_time, label = 'Bootstraps') {
    total <- as.integer(total)
    done <- as.integer(done)
    width <- 24L
    frac <- if (total > 0L) max(0, min(1, done / total)) else 1
    filled <- as.integer(floor(width * frac))
    empty <- width - filled
    elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = 'secs'))
    eta_sec <- if (done > 0L) elapsed_sec / done * (total - done) else NA_real_
    green_on <- '\033[32m'
    color_off <- '\033[39m'
    bar <- paste0(green_on, strrep('=', filled), color_off, strrep('-', empty))
    pct <- sprintf('%3d%%', round(frac * 100))
    cat(sprintf('%s [%s] %s %d/%d | elapsed %s | ETA %s\n',
                label,
                bar,
                pct,
                done,
                total,
                format_progress_time(elapsed_sec),
                format_progress_time(eta_sec)))
  }
  if (isTRUE(parallel)) {
    if (!requireNamespace('future', quietly = TRUE) || !requireNamespace('future.apply', quietly = TRUE)) {
      stop("Parallel bootstrap requires the 'future' and 'future.apply' packages.")
    }
    future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(future::sequential), add = TRUE)

    remaining_ids <- if (next_boot <= n_boot) seq.int(next_boot, n_boot) else integer()
    batches <- split_parallel_bootstrap_ids(
      remaining_ids = remaining_ids,
      n_workers = n_workers,
      batch_size = batch_size
    )
    for (i in seq_along(batches)) {
      batch_ids <- batches[[i]]
      if (length(batch_ids) == 0L) next
      batch_res <- future.apply::future_lapply(
        batch_ids,
        run_single_bootstrap,
        dt = dt,
        config = config,
        ids = ids,
        level_vars = level_vars,
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
        keep_boot_objects = keep_boot_objects,
        keep_boot_curves = keep_boot_curves,
        dt_threads = 1L,
        future.seed = TRUE,
        future.packages = c('data.table', 'splines', 'stats', 'tvcQGcomp', 'speedglm', 'nnet', 'MASS', 'qgcomp'),
        future.globals = list(
          run_single_bootstrap = run_single_bootstrap
        ),
        future.chunk.size = 1L,
        future.scheduling = Inf
      )
      for (j in seq_along(batch_ids)) {
        results[[batch_ids[j]]] <- batch_res[[j]]$result
        if (keep_boot_objects) audit_objects[[batch_ids[j]]] <- batch_res[[j]]$object
        if (keep_boot_curves) curve_results[[batch_ids[j]]] <- batch_res[[j]]$curves
      }
      rm(batch_res)
      gc(verbose = FALSE)
      save_checkpoint(max(batch_ids))
      batch_status <- unlist(lapply(results[batch_ids], function(x) x$status[[1]]), use.names = FALSE)
      err_idx <- which(grepl("^error:", batch_status))
      if (length(err_idx) > 0L && isTRUE(stop_on_error)) {
        first_err_boot <- batch_ids[[err_idx[[1L]]]]
        stop("Bootstrap ", first_err_boot, " failed: ", sub("^error:\\s*", "", batch_status[[err_idx[[1L]]]]))
      }
      if (verbose) {
        print_progress_bar(max(batch_ids), n_boot, start_time)
      }
    }
  } else {
    if (next_boot > n_boot && isTRUE(verbose)) {
      print_progress_bar(n_boot, n_boot, start_time)
    }
    for (b in seq.int(next_boot, n_boot)) {
      run_res <- run_single_bootstrap(
        b = b,
        dt = dt,
        config = config,
        ids = ids,
        level_vars = level_vars,
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
        keep_boot_objects = keep_boot_objects,
        keep_boot_curves = keep_boot_curves
      )
      results[[b]] <- run_res$result
      if (keep_boot_objects) audit_objects[[b]] <- run_res$object
      if (keep_boot_curves) curve_results[[b]] <- run_res$curves
      rm(run_res)
      gc(verbose = FALSE)
      save_checkpoint(b)
      status_b <- results[[b]]$status[[1]]
      if (grepl("^error:", status_b) && isTRUE(stop_on_error)) {
        stop("Bootstrap ", b, " failed: ", sub("^error:\\s*", "", status_b))
      }
      if (verbose) {
        print_progress_bar(b, n_boot, start_time)
      }
    }
  }

  result_dt <- data.table::rbindlist(results, fill = TRUE)
  if (!isTRUE(stop_on_error) && "status" %in% names(result_dt) && any(grepl("^error:", result_dt$status))) {
    warning(sum(grepl("^error:", result_dt$status)), " bootstrap replicate(s) failed. Inspect the 'status' column.")
  }
  if (!is.null(checkpoint_file)) {
    if (keep_boot_objects || keep_boot_curves) {
      payload <- list(results = result_dt, config = config)
      if (keep_boot_objects) payload$audit_objects <- audit_objects
      if (keep_boot_curves) payload$curve_results <- curve_results
      saveRDS(payload, checkpoint_file)
    } else {
      saveRDS(result_dt, checkpoint_file)
    }
  }

  if (keep_boot_objects || keep_boot_curves) {
    out <- list(results = result_dt[], audit_mode = isTRUE(audit_mode), config = config)
    if (keep_boot_objects) out$audit_objects <- audit_objects
    if (keep_boot_curves) out$curve_results <- data.table::rbindlist(curve_results, fill = TRUE)
    return(out)
  }

  result_dt[]
}
summarize_tvcqgcomp_bootstrap <- function(boot_res) {
  br <- if (is.list(boot_res) && !data.table::is.data.table(boot_res)) boot_res$results else boot_res
  br <- data.table::as.data.table(data.table::copy(br))
  br <- br[status == 'success']
  coef_cols <- setdiff(names(br), c('boot', 'status'))
  if (!length(coef_cols)) stop('No successful bootstrap coefficient columns were found.')

  data.table::rbindlist(lapply(coef_cols, function(v) {
    x <- br[[v]]
    qq <- stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    data.table::data.table(
      term = v,
      estimate = mean(x, na.rm = TRUE),
      se = stats::sd(x, na.rm = TRUE),
      ci_ll_95 = qq[[1]],
      ci_ul_95 = qq[[2]],
      n_success = sum(is.finite(x))
    )
  }), fill = TRUE)
}
diagnose_tvcqgcomp_run <- function(run_obj, config = NULL, tolerance = 1e-12) {
  if (is.null(config)) {
    config <- run_obj$config
  }
  if (is.null(config)) {
    stop('Please supply config or a run object that contains config.')
  }

  check_prob_cols <- function(dt, id_col, extra_group_cols = NULL) {
    group_cols <- c(id_col, intersect(extra_group_cols, names(dt)))
    out <- list(
      n_missing_hazard = if ('hazard' %in% names(dt)) sum(!is.finite(dt$hazard) | is.na(dt$hazard)) else NA_integer_,
      n_missing_cumrisk = if ('cumrisk' %in% names(dt)) sum(!is.finite(dt$cumrisk) | is.na(dt$cumrisk)) else NA_integer_,
      hazard_in_unit_interval = if ('hazard' %in% names(dt)) all(dt$hazard >= -tolerance & dt$hazard <= 1 + tolerance, na.rm = TRUE) else NA,
      cumrisk_in_unit_interval = if ('cumrisk' %in% names(dt)) all(dt$cumrisk >= -tolerance & dt$cumrisk <= 1 + tolerance, na.rm = TRUE) else NA,
      cumrisk_monotone = if ('cumrisk' %in% names(dt)) {
        dt_check <- data.table::copy(dt)
        data.table::setorderv(dt_check, c(group_cols, config$time_in))
        all(dt_check[, all(diff(cumrisk) >= -tolerance), by = group_cols]$V1)
      } else {
        NA
      }
    )
    data.table::as.data.table(out)
  }

  mc <- data.table::as.data.table(data.table::copy(run_obj$mc_data))
  nat <- if (is.null(run_obj$natural)) NULL else data.table::as.data.table(data.table::copy(run_obj$natural))
  intdat <- data.table::as.data.table(data.table::copy(run_obj$intervention_data))

  uid_col <- if ('uid' %in% names(mc)) 'uid' else config$id
  mc_counts <- mc[, .N, by = uid_col]
  mc_summary <- data.table::data.table(
    n_ids = data.table::uniqueN(mc[[uid_col]]),
    min_rows_per_id = min(mc_counts$N),
    max_rows_per_id = max(mc_counts$N),
    duplicate_id_time_rows = mc[, .N, by = c(uid_col, config$time_in)][N > 1L, .N],
    baseline_fixed_covariates = all(mc[, all(vapply(.SD, data.table::uniqueN, integer(1)) == 1L), by = uid_col, .SDcols = intersect(config$baseline_covariates, names(mc))]$V1)
  )

  natural_summary <- check_prob_cols(nat, if ('uid' %in% names(nat)) 'uid' else config$id)
  intervention_summary <- check_prob_cols(intdat, 'uid', extra_group_cols = 'scenario')

  scenario_exposure_checks <- data.table::rbindlist(lapply(names(run_obj$interventions), function(scn) {
    dt <- data.table::as.data.table(data.table::copy(run_obj$interventions[[scn]]))
    vals <- dt[, lapply(.SD, data.table::uniqueN), .SDcols = intersect(config$exposures, names(dt))]
    data.table::data.table(
      scenario = scn,
      all_exposures_constant = all(unlist(vals) == 1L),
      joint_effect_unique = if ('joint_effect' %in% names(dt)) data.table::uniqueN(dt$joint_effect) else NA_integer_
    )
  }), fill = TRUE)

  final_time <- max(intdat[[config$time_in]], na.rm = TRUE)
  final_risk_table <- intdat[get(config$time_in) == final_time,
    .(final_risk = mean(cumrisk, na.rm = TRUE), joint_effect = mean(joint_effect, na.rm = TRUE)),
    by = scenario][order(joint_effect)]
  final_risk_table[, risk_rank := data.table::frank(final_risk, ties.method = 'dense')]

  ordered_by_joint <- TRUE
  if (nrow(final_risk_table) > 1L) {
    ordered_by_joint <- all(diff(final_risk_table$final_risk) <= tolerance) || all(diff(final_risk_table$final_risk) >= -tolerance)
  }

  summary_flags <- data.table::data.table(
    mc_ok = mc_summary$duplicate_id_time_rows[[1]] == 0 && mc_summary$min_rows_per_id[[1]] == mc_summary$max_rows_per_id[[1]] && isTRUE(mc_summary$baseline_fixed_covariates[[1]]),
    natural_ok = if (is.null(nat)) NA else natural_summary$n_missing_hazard[[1]] == 0 && natural_summary$n_missing_cumrisk[[1]] == 0 && isTRUE(natural_summary$hazard_in_unit_interval[[1]]) && isTRUE(natural_summary$cumrisk_in_unit_interval[[1]]) && isTRUE(natural_summary$cumrisk_monotone[[1]]),
    intervention_ok = intervention_summary$n_missing_hazard[[1]] == 0 && intervention_summary$n_missing_cumrisk[[1]] == 0 && isTRUE(intervention_summary$hazard_in_unit_interval[[1]]) && isTRUE(intervention_summary$cumrisk_in_unit_interval[[1]]) && isTRUE(intervention_summary$cumrisk_monotone[[1]]) && all(scenario_exposure_checks$all_exposures_constant),
    scenario_order_consistent = ordered_by_joint
  )

  list(
    summary = summary_flags,
    mc = mc_summary,
    natural = natural_summary,
    intervention = intervention_summary,
    scenario_exposure_checks = scenario_exposure_checks,
    final_risk_table = final_risk_table
  )
}
diagnose_tvcqgcomp_drift <- function(observed_data = NULL,
                                     natural = NULL,
                                     run_obj = NULL,
                                     config = NULL,
                                     variables = NULL,
                                     include = c("both", "categorical", "continuous"),
                                     continuous_probs = c(0.1, 0.5, 0.9)) {
  include <- match.arg(include)

  if (!is.null(run_obj)) {
    if (is.null(config) && !is.null(run_obj$config)) config <- run_obj$config
    if (is.null(observed_data) && !is.null(run_obj$observed_data)) observed_data <- run_obj$observed_data
    if (is.null(observed_data) && !is.null(run_obj$natural_observed)) observed_data <- run_obj$natural_observed
    if (is.null(natural) && !is.null(run_obj$natural)) natural <- run_obj$natural
  }

  if (is.null(config)) {
    stop('Please supply config or a run object that contains config.')
  }
  if (is.null(observed_data) || is.null(natural)) {
    stop('Please supply observed_data and natural, or pass a run_tvcqgcomp result via run_obj.')
  }

  obs <- data.table::as.data.table(data.table::copy(observed_data))
  nat <- data.table::as.data.table(data.table::copy(natural))
  time_col <- config$time_in
  if (!(time_col %in% names(obs)) || !(time_col %in% names(nat))) {
    stop('Time column must be present in both observed_data and natural: ', time_col)
  }

  if (is.null(variables)) {
    variables <- config$tvc %||% config$time_varying_covariates
  }
  variables <- intersect(variables, intersect(names(obs), names(nat)))
  if (!length(variables)) {
    stop('No overlapping variables were found between observed_data and natural.')
  }

  probs <- sort(unique(as.numeric(continuous_probs)))
  probs <- probs[is.finite(probs) & probs >= 0 & probs <= 1]
  if (!length(probs)) {
    stop('continuous_probs must contain at least one probability between 0 and 1.')
  }

  infer_type <- function(v) {
    cfg_type <- config$covtypes[v]
    if (is.null(cfg_type) || is.na(cfg_type)) cfg_type <- config$exposure_types[v]
    cfg_type <- tolower(cfg_type %||% NA_character_)
    if (!is.na(cfg_type) && cfg_type %in% c('categorical', 'multinomial', 'binary', 'binomial', 'ordinal', 'ordinal_logit', 'ordinal_polr', 'ordered_logit', 'polr')) return('categorical')
    if (!is.na(cfg_type) && cfg_type %in% c('normal', 'continuous', 'bounded normal', 'bounded_normal', 'bounded normal snap', 'bounded_normal_snap', 'bounded normal ordinal', 'bounded_normal_ordinal')) return('continuous')
    x <- obs[[v]]
    if (is.factor(x) || is.character(x) || is.logical(x)) return('categorical')
    if (is.integer(x) && data.table::uniqueN(stats::na.omit(x)) <= 10L) return('categorical')
    'continuous'
  }

  variable_types <- data.table::data.table(
    variable = variables,
    drift_type = vapply(variables, infer_type, character(1))
  )

  categorical_vars <- variable_types[drift_type == 'categorical', variable]
  continuous_vars <- variable_types[drift_type == 'continuous', variable]
  if (identical(include, 'categorical')) continuous_vars <- character(0)
  if (identical(include, 'continuous')) categorical_vars <- character(0)

  categorical_detail <- data.table::data.table()
  categorical_summary <- data.table::data.table()
  if (length(categorical_vars)) {
    categorical_detail <- data.table::rbindlist(lapply(categorical_vars, function(v) {
      obs_counts <- obs[, .N, by = c(time_col, v)]
      data.table::setnames(obs_counts, v, 'value')
      obs_counts[, p_observed := N / sum(N), by = c(time_col)]
      nat_counts <- nat[, .N, by = c(time_col, v)]
      data.table::setnames(nat_counts, v, 'value')
      nat_counts[, p_natural := N / sum(N), by = c(time_col)]
      out <- merge(
        obs_counts[, .(time_value = get(time_col), value, p_observed)],
        nat_counts[, .(time_value = get(time_col), value, p_natural)],
        by = c('time_value', 'value'),
        all = TRUE
      )
      out[is.na(p_observed), p_observed := 0]
      out[is.na(p_natural), p_natural := 0]
      out[, `:=`(
        variable = v,
        signed_diff = p_natural - p_observed,
        abs_diff = abs(p_natural - p_observed)
      )]
      data.table::setcolorder(out, c('variable', 'time_value', 'value', 'p_observed', 'p_natural', 'signed_diff', 'abs_diff'))
      out[order(time_value, value)]
    }), fill = TRUE)

    categorical_summary <- categorical_detail[, .(
      n_time = data.table::uniqueN(time_value),
      n_levels = data.table::uniqueN(value),
      max_abs_diff = max(abs_diff, na.rm = TRUE),
      p95_abs_diff = stats::quantile(abs_diff, probs = 0.95, na.rm = TRUE, names = FALSE),
      mean_abs_diff = mean(abs_diff, na.rm = TRUE)
    ), by = variable][order(-max_abs_diff, variable)]
  }

  continuous_detail <- data.table::data.table()
  continuous_quantiles <- data.table::data.table()
  continuous_summary <- data.table::data.table()
  if (length(continuous_vars)) {
    fmt_prob <- function(p) paste0('q', gsub('\\.', '_', format(round(p * 100, 4), trim = TRUE, scientific = FALSE)))

    continuous_detail <- data.table::rbindlist(lapply(continuous_vars, function(v) {
      obs_stats <- obs[, .(
        n_observed = sum(!is.na(get(v))),
        mean_observed = mean(get(v), na.rm = TRUE),
        sd_observed = stats::sd(get(v), na.rm = TRUE)
      ), by = c(time_col)]
      nat_stats <- nat[, .(
        n_natural = sum(!is.na(get(v))),
        mean_natural = mean(get(v), na.rm = TRUE),
        sd_natural = stats::sd(get(v), na.rm = TRUE)
      ), by = c(time_col)]
      out <- merge(obs_stats, nat_stats, by = time_col, all = TRUE)
      out[, `:=`(
        variable = v,
        mean_diff = mean_natural - mean_observed,
        abs_mean_diff = abs(mean_natural - mean_observed),
        sd_diff = sd_natural - sd_observed,
        abs_sd_diff = abs(sd_natural - sd_observed)
      )]
      data.table::setnames(out, time_col, 'time_value')
      data.table::setcolorder(out, c('variable', 'time_value', 'n_observed', 'n_natural', 'mean_observed', 'mean_natural', 'mean_diff', 'abs_mean_diff', 'sd_observed', 'sd_natural', 'sd_diff', 'abs_sd_diff'))
      out[order(time_value)]
    }), fill = TRUE)

    continuous_quantiles <- data.table::rbindlist(lapply(continuous_vars, function(v) {
      obs_q <- obs[, as.list(stats::quantile(get(v), probs = probs, na.rm = TRUE, names = FALSE)), by = c(time_col)]
      nat_q <- nat[, as.list(stats::quantile(get(v), probs = probs, na.rm = TRUE, names = FALSE)), by = c(time_col)]
      qnames <- vapply(probs, fmt_prob, character(1))
      data.table::setnames(obs_q, old = names(obs_q)[-1], new = qnames)
      data.table::setnames(nat_q, old = names(nat_q)[-1], new = qnames)
      out <- merge(obs_q, nat_q, by = time_col, all = TRUE, suffixes = c('_observed', '_natural'))
      out_long <- data.table::rbindlist(lapply(seq_along(probs), function(i) {
        qn <- qnames[[i]]
        obs_col <- paste0(qn, '_observed')
        nat_col <- paste0(qn, '_natural')
        data.table::data.table(
          variable = v,
          time_value = out[[time_col]],
          quantile = probs[[i]],
          value_observed = out[[obs_col]],
          value_natural = out[[nat_col]],
          diff = out[[nat_col]] - out[[obs_col]],
          abs_diff = abs(out[[nat_col]] - out[[obs_col]])
        )
      }), fill = TRUE)
      out_long[order(time_value, quantile)]
    }), fill = TRUE)

    quant_summary <- if (nrow(continuous_quantiles)) {
      continuous_quantiles[, .(
        max_abs_quantile_diff = max(abs_diff, na.rm = TRUE),
        mean_abs_quantile_diff = mean(abs_diff, na.rm = TRUE)
      ), by = variable]
    } else {
      data.table::data.table(variable = continuous_vars, max_abs_quantile_diff = NA_real_, mean_abs_quantile_diff = NA_real_)
    }

    continuous_summary <- merge(
      continuous_detail[, .(
        n_time = data.table::uniqueN(time_value),
        max_abs_mean_diff = max(abs_mean_diff, na.rm = TRUE),
        mean_abs_mean_diff = mean(abs_mean_diff, na.rm = TRUE),
        max_abs_sd_diff = max(abs_sd_diff, na.rm = TRUE),
        mean_abs_sd_diff = mean(abs_sd_diff, na.rm = TRUE)
      ), by = variable],
      quant_summary,
      by = 'variable',
      all = TRUE
    )[order(-max_abs_mean_diff, variable)]
  }

  list(
    settings = list(
      time_col = time_col,
      variables = variables,
      include = include,
      continuous_probs = probs
    ),
    variable_types = variable_types,
    categorical_detail = categorical_detail,
    categorical_summary = categorical_summary,
    continuous_detail = continuous_detail,
    continuous_quantiles = continuous_quantiles,
    continuous_summary = continuous_summary
  )
}

plot_tvcqgcomp_drift <- function(drift_obj = NULL,
                                 observed_data = NULL,
                                 natural = NULL,
                                 run_obj = NULL,
                                 config = NULL,
                                 variables = NULL,
                                 include = c("both", "categorical", "continuous"),
                                 continuous_stat = c("mean", "sd", "quantile"),
                                 quantile = 0.5,
                                 categorical_stat = c("mean_score", "level_proportion"),
                                 categorical_level = NULL,
                                 ncol = NULL,
                                 observed_col = "black",
                                 natural_col = "grey70",
                                 observed_pch = 16,
                                 natural_lwd = 2,
                                 xlab = NULL,
                                 ylab = NULL,
                                 main = NULL,
                                 legend = TRUE,
                                 legend_pos = "bottom",
                                 ...) {
  include <- match.arg(include)
  continuous_stat <- match.arg(continuous_stat)
  categorical_stat <- match.arg(categorical_stat)

  if (is.null(drift_obj)) {
    drift_obj <- diagnose_tvcqgcomp_drift(
      observed_data = observed_data,
      natural = natural,
      run_obj = run_obj,
      config = config,
      variables = variables,
      include = include,
      continuous_probs = quantile
    )
  }

  panels <- list()

  if (include %in% c("both", "continuous") && nrow(drift_obj$continuous_detail)) {
    if (identical(continuous_stat, "quantile")) {
      q_target <- quantile
      qdat <- drift_obj$continuous_quantiles[abs(quantile - q_target) < 1e-10]
      if (!nrow(qdat)) {
        stop("Requested quantile was not found in the drift object. Re-run diagnose_tvcqgcomp_drift() with continuous_probs including ", q_target, ".")
      }
      cont_panels <- lapply(split(qdat, by = "variable", keep.by = FALSE), function(dat) {
        dat <- data.table::as.data.table(dat)
        dat[, observed_value := value_observed]
        dat[, natural_value := value_natural]
        data.table::setorderv(dat, "time_value")
        dat[, .(time_value, observed_value, natural_value)]
      })
    } else {
      cdat <- data.table::copy(drift_obj$continuous_detail)
      obs_col <- if (identical(continuous_stat, "mean")) "mean_observed" else "sd_observed"
      nat_col <- if (identical(continuous_stat, "mean")) "mean_natural" else "sd_natural"
      cont_panels <- lapply(split(cdat, by = "variable", keep.by = FALSE), function(dat) {
        dat <- data.table::as.data.table(dat)
        dat[, observed_value := get(obs_col)]
        dat[, natural_value := get(nat_col)]
        data.table::setorderv(dat, "time_value")
        dat[, .(time_value, observed_value, natural_value)]
      })
    }
    panels <- c(panels, cont_panels)
  }

  if (include %in% c("both", "categorical") && nrow(drift_obj$categorical_detail)) {
    cat_panels <- lapply(split(drift_obj$categorical_detail, by = "variable", keep.by = FALSE), function(dat) {
      dat <- data.table::as.data.table(dat)
      data.table::setorderv(dat, c("time_value", "value"))
      if (identical(categorical_stat, "level_proportion")) {
        if (is.null(categorical_level)) {
          top_level <- dat[, .(mean_observed = mean(p_observed, na.rm = TRUE)), by = value][order(-mean_observed)][1L, value]
        } else {
          top_level <- categorical_level
        }
        out <- dat[value == top_level, .(
          time_value,
          observed_value = p_observed,
          natural_value = p_natural
        )]
        attr(out, "subtitle") <- paste0("Level: ", top_level)
        return(out)
      }

      suppressWarnings(value_num <- as.numeric(as.character(dat$value)))
      if (anyNA(value_num)) {
        levels_in_order <- unique(as.character(dat$value))
        value_num <- match(as.character(dat$value), levels_in_order)
      }
      dat[, value_num := value_num]
      out <- dat[, .(
        observed_value = sum(value_num * p_observed, na.rm = TRUE),
        natural_value = sum(value_num * p_natural, na.rm = TRUE)
      ), by = time_value]
      out
    })
    panels <- c(panels, cat_panels)
  }

  if (!length(panels)) {
    stop("No drift panels available for plotting with the requested settings.")
  }

  panel_names <- names(panels)
  n_panels <- length(panels)
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_panels))
  }
  nrow <- ceiling(n_panels / ncol)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(nrow, ncol), mar = c(3.2, 3.4, 3.6, 1.2), oma = c(0, 0, if (!is.null(main)) 2 else 0, if (legend) 2 else 0))

  for (i in seq_along(panels)) {
    dat <- data.table::as.data.table(panels[[i]])
    data.table::setorderv(dat, "time_value")
    ylim <- range(c(dat$observed_value, dat$natural_value), na.rm = TRUE)
    if (!all(is.finite(ylim))) ylim <- c(0, 1)
    if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
    subtitle <- attr(panels[[i]], "subtitle", exact = TRUE)
    panel_title <- if (!is.null(subtitle)) paste0(panel_names[[i]], "\n", subtitle) else panel_names[[i]]
    graphics::plot(dat$time_value, dat$natural_value,
      type = "l",
      col = natural_col,
      lwd = natural_lwd,
      xlab = if (is.null(xlab)) "" else xlab,
      ylab = if (is.null(ylab)) "" else ylab,
      main = panel_title,
      ylim = ylim,
      ...
    )
    graphics::points(dat$time_value, dat$observed_value, pch = observed_pch, col = observed_col)
  }

  if (!is.null(main)) {
    graphics::mtext(main, side = 3, outer = TRUE, line = 0.5, cex = 1.1)
  }
  if (legend) {
    graphics::par(xpd = NA)
    graphics::legend(
      legend_pos,
      inset = -0.02,
      legend = c("Observed", "Natural Course"),
      col = c(observed_col, natural_col),
      pch = c(observed_pch, NA),
      lwd = c(NA, natural_lwd),
      horiz = TRUE,
      bty = "n"
    )
  }

  invisible(drift_obj)
}

plot_tvc_drift <- function(observed_data,
                           run_obj,
                                      natural_source = NULL,
                                      time_var = NULL,
                                      exposure_vars = NULL,
                                      categorical_tvc_vars = NULL,
                                      continuous_tvc_vars = NULL,
                                      bounded_normal_vars = NULL,
                                      bounded_normal_summary = c("median", "mean"),
                                      observed_label = "Observed",
                                      natural_label = "Natural Course",
                                      ncol_main = 3,
                                      x_breaks = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_tvc_drift() requires the 'ggplot2' package. Install it with install.packages('ggplot2').")
  }
  bounded_normal_summary <- match.arg(bounded_normal_summary)

  config   <- run_obj$config
  if (is.null(time_var))      time_var      <- config$time_in
  if (is.null(exposure_vars)) exposure_vars <- config$exposures

  # Resolve natural course source
  if (is.null(natural_source)) {
    if (!is.null(run_obj$calibration$modeled_exposome))   natural_source <- "modeled_exposome"
    else if (!is.null(run_obj$calibration$observed_exposome)) natural_source <- "observed_exposome"
    else if (!is.null(run_obj$calibration$calibration))   natural_source <- "calibration"
    else natural_source <- "modeled_exposome"
  }
  nat_data <- run_obj$calibration[[natural_source]]
  if (is.null(nat_data)) stop("Could not find natural course data for source: ", natural_source)

  obs <- data.table::as.data.table(data.table::copy(observed_data))
  nat <- data.table::as.data.table(data.table::copy(nat_data))

  # If the config used quantized exposures, bring observed exposures onto the same
  # quantile-index scale so that observed vs natural-course means are comparable.
  exposure_scale <- tolower(config$exposome %||% config$exposure_scale %||% "quantized")
  if (identical(exposure_scale, "quantized") && length(exposure_vars) > 0L) {
    obs_q <- tryCatch(
      quantize_exposures(obs, config),
      error = function(e) NULL
    )
    if (!is.null(obs_q)) {
      for (v in intersect(exposure_vars, names(obs_q))) {
        obs[, (v) := obs_q[[v]]]
      }
    }
  }

  # Auto-detect variable types from config if not supplied
  if (is.null(categorical_tvc_vars) && is.null(continuous_tvc_vars)) {
    categorical_tvc_vars <- names(config$covtypes)[tolower(config$covtypes) %in%
      c("categorical", "multinomial", "binary", "binomial", "ordinal", "ordinal_logit", "ordinal_polr", "ordered_logit", "polr")]
    continuous_tvc_vars  <- names(config$covtypes)[tolower(config$covtypes) %in%
      c("normal", "continuous", "bounded normal", "bounded_normal", "bounded normal snap", "bounded_normal_snap", "bounded normal ordinal", "bounded_normal_ordinal")]
  }
  auto_bounded_exposures <- names(config$exposure_types)[tolower(config$exposure_types) %in%
    c("bounded normal", "bounded_normal", "bounded normal snap", "bounded_normal_snap", "bounded normal ordinal", "bounded_normal_ordinal")]
  auto_bounded_tvcs <- names(config$covtypes)[tolower(config$covtypes) %in%
    c("bounded normal", "bounded_normal", "bounded normal snap", "bounded_normal_snap", "bounded normal ordinal", "bounded_normal_ordinal")]
  if (is.null(bounded_normal_vars)) {
    bounded_normal_vars <- unique(c(auto_bounded_exposures, auto_bounded_tvcs))
  }
  categorical_tvc_vars <- intersect(categorical_tvc_vars %||% character(), intersect(names(obs), names(nat)))
  continuous_tvc_vars  <- intersect(continuous_tvc_vars  %||% character(), intersect(names(obs), names(nat)))
  bounded_normal_vars  <- intersect(bounded_normal_vars  %||% character(), intersect(names(obs), names(nat)))
  continuous_tvc_vars  <- setdiff(continuous_tvc_vars, bounded_normal_vars)
  exposure_vars        <- intersect(exposure_vars %||% character(), intersect(names(obs), names(nat)))
  bounded_exposure_vars <- intersect(exposure_vars, bounded_normal_vars)
  mean_exposure_vars <- setdiff(exposure_vars, bounded_normal_vars)

  # Align factor levels for categorical TVCs
  for (v in categorical_tvc_vars) {
    if (!is.factor(obs[[v]])) obs[, (v) := factor(get(v))]
    nat[, (v) := factor(as.character(get(v)), levels = levels(obs[[v]]))]
  }

  # Mean over time for continuous/exposure variables
  make_mean_long <- function(dt, vars, source_label, panel_type) {
    if (!length(vars)) return(NULL)
    x <- data.table::copy(dt[, c(time_var, vars), with = FALSE])
    for (v in vars) {
      if (is.factor(x[[v]])) x[, (v) := as.numeric(as.character(get(v)))]
    }
    out <- data.table::melt(x, id.vars = time_var, variable.name = "variable", value.name = "value")
    out <- out[, .(value = mean(value, na.rm = TRUE)), by = c(time_var, "variable")]
    out[, `:=`(source = source_label, level = NA_character_, panel_type = panel_type)]
    out
  }

  make_bounded_summary_long <- function(obs_dt, nat_dt, vars, obs_label, nat_label, panel_type, summary_stat) {
    if (!length(vars)) return(NULL)
    data.table::rbindlist(lapply(vars, function(v) {
      obs_num <- suppressWarnings(as.numeric(as.character(obs_dt[[v]])))
      nat_num <- suppressWarnings(as.numeric(as.character(nat_dt[[v]])))
      finite_obs <- obs_num[is.finite(obs_num)]
      if (!length(finite_obs)) return(NULL)
      min_level <- floor(min(finite_obs, na.rm = TRUE))
      max_level <- ceiling(max(finite_obs, na.rm = TRUE))

      obs_x <- data.table::copy(obs_dt[, c(time_var), with = FALSE])
      nat_x <- data.table::copy(nat_dt[, c(time_var), with = FALSE])
      obs_x[, value_tmp := pmin(pmax(round(obs_num), min_level), max_level)]
      nat_x[, value_tmp := pmin(pmax(round(nat_num), min_level), max_level)]

      obs_out <- if (identical(summary_stat, "mean")) {
        obs_x[, .(value = mean(value_tmp, na.rm = TRUE)), by = c(time_var)]
      } else {
        obs_x[, .(value = stats::median(value_tmp, na.rm = TRUE)), by = c(time_var)]
      }
      nat_out <- if (identical(summary_stat, "mean")) {
        nat_x[, .(value = mean(value_tmp, na.rm = TRUE)), by = c(time_var)]
      } else {
        nat_x[, .(value = stats::median(value_tmp, na.rm = TRUE)), by = c(time_var)]
      }
      obs_out[, `:=`(variable = v, source = obs_label, level = NA_character_, panel_type = panel_type)]
      nat_out[, `:=`(variable = v, source = nat_label, level = NA_character_, panel_type = panel_type)]
      data.table::rbindlist(list(obs_out, nat_out), fill = TRUE)
    }), fill = TRUE)
  }

  # Proportion of every level over time for categorical variables
  make_prop_long <- function(dt, vars, source_label) {
    if (!length(vars)) return(NULL)
    data.table::rbindlist(lapply(vars, function(v) {
      x <- data.table::copy(dt[, c(time_var, v), with = FALSE])
      if (!is.factor(x[[v]])) x[, (v) := factor(get(v))]
      levs    <- levels(x[[v]])
      counts  <- x[, .N, by = c(time_var, v)]
      data.table::setnames(counts, v, "level")
      totals  <- x[, .(denom = .N), by = time_var]
      counts  <- merge(counts, totals, by = time_var, all.x = TRUE)
      counts[, value := N / denom]
      counts[, `:=`(variable = v, source = source_label, panel_type = "categorical_tvc")]
      # Ensure every time × level combination is present (fill zeros)
      full <- data.table::CJ(tmp_time = sort(unique(x[[time_var]])), tmp_level = levs)
      data.table::setnames(full, c("tmp_time", "tmp_level"), c(time_var, "level"))
      full[, level := factor(level, levels = levs)]
      counts[, level := factor(as.character(level), levels = levs)]
      out <- merge(full,
                   counts[, c(time_var, "level", "value", "variable", "source", "panel_type"), with = FALSE],
                   by = c(time_var, "level"), all.x = TRUE)
      out[is.na(value),      value      := 0]
      out[is.na(variable),   variable   := v]
      out[is.na(source),     source     := source_label]
      out[is.na(panel_type), panel_type := "categorical_tvc"]
      out
    }), fill = TRUE)
  }

  # Survival curves: compare empirical observed survival to natural-course survival.
  make_empirical_survival_curve <- function(dt_in, cfg) {
    dt_surv <- data.table::as.data.table(data.table::copy(dt_in))
    id_var <- cfg$id
    outcome_var <- cfg$outcome
    if (!(id_var %in% names(dt_surv)) || !(time_var %in% names(dt_surv)) || !(outcome_var %in% names(dt_surv))) {
      return(NULL)
    }
    data.table::setorderv(dt_surv, c(id_var, time_var))
    dt_surv[, .prior_event := data.table::shift(cumsum(get(outcome_var) > 0), fill = 0L), by = id_var]
    risk_dt <- dt_surv[.prior_event == 0L]
    if (!nrow(risk_dt)) {
      return(NULL)
    }
    surv_dt <- risk_dt[, .(
      n_risk = data.table::uniqueN(get(id_var)),
      n_event = sum(get(outcome_var) > 0, na.rm = TRUE)
    ), by = c(time_var)]
    surv_dt[, hazard := data.table::fifelse(n_risk > 0, n_event / n_risk, NA_real_)]
    surv_dt[, survival := cumprod(1 - hazard)]
    surv_dt[, .(value = survival), by = c(time_var)]
  }

  obs_surv <- make_empirical_survival_curve(observed_data, config)
  nat_surv_traj <- summarize_risk_trajectory(nat_data, config = config)
  nat_surv <- nat_surv_traj[, .(value = mean_survival), by = c(time_var)]
  nat_surv[, `:=`(variable = "Survival", source = natural_label,  level = NA_character_, panel_type = "survival")]
  if (!is.null(obs_surv)) {
    obs_surv[, `:=`(variable = "Survival", source = observed_label, level = NA_character_, panel_type = "survival")]
    surv_cmp <- merge(
      obs_surv[, c(time_var, "value"), with = FALSE],
      nat_surv[, c(time_var, "value"), with = FALSE],
      by = time_var,
      suffixes = c("_obs", "_nat"),
      all = FALSE
    )
    identical_survival <- nrow(surv_cmp) > 0L &&
      isTRUE(all.equal(surv_cmp$value_obs, surv_cmp$value_nat, tolerance = 0))
    if (identical_survival) {
      warning(
        "Survival panel omitted because observed/calibration and natural-course survival trajectories are identical.",
        call. = FALSE
      )
      obs_surv <- NULL
      nat_surv <- NULL
    }
  } else {
    warning(
      "No empirical observed survival trajectory could be constructed; survival panel omitted.",
      call. = FALSE
    )
    nat_surv <- NULL
  }

  plot_dt <- data.table::rbindlist(list(
    obs_surv, nat_surv,
    make_mean_long(obs, mean_exposure_vars,        observed_label, "exposure"),
    make_mean_long(nat, mean_exposure_vars,        natural_label,  "exposure"),
    make_bounded_summary_long(obs, nat, bounded_exposure_vars, observed_label, natural_label, "bounded_normal", bounded_normal_summary),
    make_mean_long(obs, continuous_tvc_vars,  observed_label, "continuous_tvc"),
    make_mean_long(nat, continuous_tvc_vars,  natural_label,  "continuous_tvc"),
    make_bounded_summary_long(obs, nat, setdiff(bounded_normal_vars, bounded_exposure_vars), observed_label, natural_label, "bounded_normal", bounded_normal_summary),
    make_prop_long(obs, categorical_tvc_vars, observed_label),
    make_prop_long(nat, categorical_tvc_vars, natural_label)
  ), fill = TRUE)

  plot_dt[, source := factor(source, levels = c(observed_label, natural_label))]

  if (is.null(x_breaks)) {
    times    <- sort(unique(plot_dt[[time_var]]))
    x_breaks <- if (length(times) <= 15) times else pretty(times, n = 10)
  }

  pal <- c("#D55E00", "#0072B2")
  names(pal) <- c(observed_label, natural_label)
  shape_map <- c(1, 17)
  names(shape_map) <- c(observed_label, natural_label)
  size_map <- c(2.9, 2.3)
  names(size_map) <- c(observed_label, natural_label)

  # Main plot: survival + exposures + continuous TVCs
  main_dt <- plot_dt[panel_type %in% c("survival", "exposure", "continuous_tvc", "bounded_normal")]

  # Dummy data to force the Survival panel y-axis to [0, 1] while other panels stay free
  t0 <- min(main_dt[[time_var]], na.rm = TRUE)
  surv_bounds <- data.table::data.table(
    value    = c(0, 1),
    variable = "Survival",
    source   = observed_label
  )
  surv_bounds[, (time_var) := t0]

  bounded_bounds <- NULL
  bounded_main_vars <- intersect(unique(main_dt[panel_type == "bounded_normal", variable]), names(obs))
  if (length(bounded_main_vars) > 0L) {
    bounded_bounds <- data.table::rbindlist(lapply(bounded_main_vars, function(v) {
      vals <- suppressWarnings(as.numeric(as.character(obs[[v]])))
      vals <- vals[is.finite(vals)]
      if (!length(vals)) return(NULL)
      out <- data.table::data.table(
        value = c(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)),
        variable = v,
        source = observed_label
      )
      out[, (time_var) := t0]
      out
    }), fill = TRUE)
  }

  axis_bounds <- data.table::rbindlist(list(surv_bounds, bounded_bounds), fill = TRUE)

  main_dt_plot <- data.table::copy(main_dt)
  main_dt_plot[, x_plot := get(time_var)]

  plot_main <- ggplot2::ggplot(main_dt_plot,
      ggplot2::aes(x = x_plot, y = value, color = source, shape = source, size = source)) +
    ggplot2::geom_blank(data = axis_bounds,
      ggplot2::aes(x = .data[[time_var]], y = value)) +
    ggplot2::geom_point(alpha = 0.95, stroke = 1.1) +
    ggplot2::facet_wrap(~ variable, scales = "free_y", ncol = ncol_main) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_shape_manual(values = shape_map) +
    ggplot2::scale_size_manual(values = size_map, guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
    ggplot2::labs(
      x = "Follow-up time", y = NULL, color = NULL, shape = NULL,
      title = paste("Observed vs", natural_label,
                    "\u2014 Survival, Exposures & Continuous TVCs")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor    = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.title          = ggplot2::element_text(size = ggplot2::rel(1.2), hjust = 0),
      axis.title.x        = ggplot2::element_text(size = ggplot2::rel(0.95)),
      axis.text           = ggplot2::element_text(size = ggplot2::rel(0.9)),
      strip.text          = ggplot2::element_text(size = ggplot2::rel(0.95)),
      legend.text         = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.position     = "bottom")

  # Categorical plot: all levels per TVC variable
  cat_dt     <- plot_dt[panel_type == "categorical_tvc"]
  plot_cat   <- NULL
  if (nrow(cat_dt) > 0L) {
    # Detect whether all categorical vars have the same set of levels
    level_sets    <- lapply(categorical_tvc_vars, function(v) {
      sort(unique(as.character(cat_dt[variable == v, level])))
    })
    levels_equal  <- length(unique(level_sets)) == 1L
    n_vars_cat    <- length(categorical_tvc_vars)

    if (levels_equal) {
      # All variables share the same levels → clean aligned grid
      facet_spec <- ggplot2::facet_grid(level ~ variable, scales = "free_y")
    } else {
      # Variables have different numbers of levels → facet_wrap avoids empty panels
      # Create a combined label ordered so level varies slowest (rows) and variable fastest (cols)
      all_levs <- sort(unique(as.character(cat_dt$level)))
      cat_dt[, .panel := factor(
        paste0(variable, "\nLevel: ", level),
        levels = as.vector(outer(all_levs, categorical_tvc_vars,
                                 function(l, v) paste0(v, "\nLevel: ", l)))
      )]
      facet_spec <- ggplot2::facet_wrap(
        ~ .panel, scales = "free_y", ncol = n_vars_cat, drop = TRUE)
    }

    cat_dt_plot <- data.table::copy(cat_dt)
    cat_dt_plot[, x_plot := get(time_var)]

    base_gg <- ggplot2::ggplot(cat_dt_plot,
        ggplot2::aes(x = x_plot, y = value, color = source, shape = source, size = source)) +
      ggplot2::geom_point(alpha = 0.95, stroke = 1.1) +
      facet_spec +
      ggplot2::scale_color_manual(values = pal) +
      ggplot2::scale_shape_manual(values = shape_map) +
      ggplot2::scale_size_manual(values = size_map, guide = "none") +
      ggplot2::scale_x_continuous(
        breaks = x_breaks,
        guide = ggplot2::guide_axis(check.overlap = TRUE)
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
      ggplot2::labs(
        x = "Follow-up time", y = "Proportion", color = NULL, shape = NULL,
        title = paste("Observed vs", natural_label,
                      "\u2014 Category Proportions for TVCs")) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor    = ggplot2::element_blank(),
        plot.title.position = "plot",
        plot.title          = ggplot2::element_text(size = ggplot2::rel(1.2), hjust = 0),
        axis.title.x        = ggplot2::element_text(size = ggplot2::rel(0.95)),
        axis.title.y        = ggplot2::element_text(size = ggplot2::rel(0.95)),
        axis.text           = ggplot2::element_text(size = ggplot2::rel(0.9)),
        strip.text          = ggplot2::element_text(size = ggplot2::rel(0.95)),
        legend.text         = ggplot2::element_text(size = ggplot2::rel(0.9)),
        legend.position     = "bottom")

    plot_cat <- base_gg
  }

  list(data = plot_dt, plot_main = plot_main, plot_categorical = plot_cat)
}

diagnose_positivity_tvcqgcomp <- function(observed_data = NULL,
                                          config,
                                          natural = NULL,
                                          intervention_data = NULL,
                                          run_obj = NULL,
                                          model_fit = NULL,
                                          scenario_col = 'scenario') {
  if (!is.null(run_obj)) {
    if (is.null(observed_data) && !is.null(run_obj$observed_data)) observed_data <- run_obj$observed_data
    if (is.null(observed_data) && !is.null(run_obj$natural_observed)) observed_data <- run_obj$natural_observed
    if (is.null(natural) && !is.null(run_obj$natural)) natural <- run_obj$natural
    if (is.null(intervention_data) && !is.null(run_obj$intervention_data)) intervention_data <- run_obj$intervention_data
    if (is.null(model_fit) && !is.null(run_obj$model_fit)) model_fit <- run_obj$model_fit
  }

  if (is.null(observed_data) || is.null(natural) || is.null(intervention_data)) {
    stop('Please supply observed_data, natural, and intervention_data, or pass a run_tvcqgcomp result via run_obj.')
  }

  obs <- data.table::as.data.table(observed_data)
  nat <- data.table::as.data.table(natural)
  int <- data.table::as.data.table(intervention_data)

  exposure_scale <- tolower(config$exposome %||% config$exposure_scale %||% "quantized")
  if (identical(exposure_scale, "quantized")) {
    q_levels <- config$q_levels %||% config$q %||% 4L
    expected_levels <- 0:(as.integer(q_levels) - 1L)
    is_already_quantized <- all(vapply(intersect(config$exposures, names(obs)), function(v) {
      x <- obs[[v]]
      x <- x[is.finite(x) & !is.na(x)]
      if (!length(x)) return(TRUE)
      all(abs(x - round(x)) < 1e-8) && all(round(x) %in% expected_levels)
    }, logical(1)))

    if (!is_already_quantized) {
      qprep <- prepare_exposure_scale_data(
        data = obs,
        config = config,
        exposure_scale = "quantized",
        q = q_levels
      )
      obs <- data.table::as.data.table(qprep$data)
    }
  }

  if (!(scenario_col %in% names(int))) {
    int[, (scenario_col) := 'intervention']
  }

  scenarios <- unique(int[[scenario_col]])

  overall_support <- data.table::rbindlist(lapply(scenarios, function(scn) {
    scn_dt <- int[get(scenario_col) == scn]
    data.table::rbindlist(lapply(config$exposures, function(v) {
      obs_min <- min(obs[[v]], na.rm = TRUE)
      obs_max <- max(obs[[v]], na.rm = TRUE)
      nat_vals <- nat[[v]]
      int_vals <- scn_dt[[v]]
      data.table::data.table(
        scenario = scn,
        exposure = v,
        obs_min = obs_min,
        obs_max = obs_max,
        natural_min = min(nat_vals, na.rm = TRUE),
        natural_max = max(nat_vals, na.rm = TRUE),
        intervention_min = min(int_vals, na.rm = TRUE),
        intervention_max = max(int_vals, na.rm = TRUE),
        prop_natural_outside = mean(nat_vals < obs_min | nat_vals > obs_max, na.rm = TRUE),
        prop_intervention_outside = mean(int_vals < obs_min | int_vals > obs_max, na.rm = TRUE)
      )
    }), fill = TRUE)
  }), fill = TRUE)

  time_support <- data.table::rbindlist(lapply(scenarios, function(scn) {
    scn_dt <- int[get(scenario_col) == scn]
    data.table::rbindlist(lapply(config$exposures, function(v) {
      obs_time <- obs[, .(obs_min = min(get(v), na.rm = TRUE), obs_max = max(get(v), na.rm = TRUE)), by = c(config$time_in)]
      nat_time <- nat[, .(natural_min = min(get(v), na.rm = TRUE), natural_max = max(get(v), na.rm = TRUE)), by = c(config$time_in)]
      int_time <- scn_dt[, .(intervention_min = min(get(v), na.rm = TRUE), intervention_max = max(get(v), na.rm = TRUE)), by = c(config$time_in)]

      nat_flag <- merge(nat[, .(time_value = get(config$time_in), value = get(v))],
                        obs_time[, .(time_value = get(config$time_in), obs_min, obs_max)],
                        by = 'time_value', all.x = TRUE)
      int_flag <- merge(scn_dt[, .(time_value = get(config$time_in), value = get(v))],
                        obs_time[, .(time_value = get(config$time_in), obs_min, obs_max)],
                        by = 'time_value', all.x = TRUE)

      nat_out <- nat_flag[, .(prop_natural_outside_time_support = mean(value < obs_min | value > obs_max, na.rm = TRUE)), by = time_value]
      int_out <- int_flag[, .(prop_intervention_outside_time_support = mean(value < obs_min | value > obs_max, na.rm = TRUE)), by = time_value]

      out <- merge(obs_time, nat_time, by = config$time_in, all = TRUE)
      out <- merge(out, int_time, by = config$time_in, all = TRUE)
      data.table::setnames(nat_out, 'time_value', config$time_in)
      data.table::setnames(int_out, 'time_value', config$time_in)
      out <- merge(out, nat_out, by = config$time_in, all = TRUE)
      out <- merge(out, int_out, by = config$time_in, all = TRUE)
      out[, `:=`(scenario = scn, exposure = v)]
      out
    }), fill = TRUE)
  }), fill = TRUE)

  mean_by_time <- data.table::rbindlist(lapply(scenarios, function(scn) {
    scn_dt <- int[get(scenario_col) == scn]
    data.table::rbindlist(lapply(config$exposures, function(v) {
      nat_mean <- nat[, .(mean_natural = mean(get(v), na.rm = TRUE)), by = c(config$time_in)]
      int_mean <- scn_dt[, .(mean_intervention = mean(get(v), na.rm = TRUE)), by = c(config$time_in)]
      out <- merge(nat_mean, int_mean, by = config$time_in, all = TRUE)
      out[, `:=`(scenario = scn, exposure = v, mean_diff = mean_intervention - mean_natural)]
      out
    }), fill = TRUE)
  }), fill = TRUE)

  joint_support_overall <- NULL
  joint_support_time <- NULL
  if (length(intersect(config$exposures, names(obs))) == length(config$exposures)) {
    obs_joint <- obs[, .N, by = c(config$exposures)]
    total_obs_n <- nrow(obs)
    obs_joint[, observed_prop := N / total_obs_n]

    scenario_joint <- unique(int[, c(scenario_col, config$exposures), with = FALSE])
    joint_support_overall <- merge(
      scenario_joint,
      obs_joint,
      by = config$exposures,
      all.x = TRUE
    )
    joint_support_overall[is.na(N), `:=`(N = 0L, observed_prop = 0)]
    data.table::setnames(joint_support_overall, "N", "observed_N")
    joint_support_overall[, observed_present := observed_N > 0L]

    obs_joint_time <- obs[, .N, by = c(config$time_in, config$exposures)]
    time_totals <- obs[, .(time_total_N = .N), by = c(config$time_in)]
    obs_joint_time <- merge(obs_joint_time, time_totals, by = config$time_in, all.x = TRUE)
    obs_joint_time[, observed_prop_time := N / time_total_N]

    scenario_joint_time <- unique(int[, c(scenario_col, config$time_in, config$exposures), with = FALSE])
    joint_support_time <- merge(
      scenario_joint_time,
      obs_joint_time,
      by = c(config$time_in, config$exposures),
      all.x = TRUE
    )
    joint_support_time[is.na(N), `:=`(N = 0L, observed_prop_time = 0)]
    data.table::setnames(joint_support_time, "N", "observed_N")
    joint_support_time[, observed_present := observed_N > 0L]
  }

  overlap_vars <- unique(intersect(c(config$time_out, config$exposures, config$baseline_covariates, config$tvc), names(obs)))
  overlap_summary <- data.table::rbindlist(lapply(scenarios, function(scn) {
    scn_dt <- int[get(scenario_col) == scn]
    if (length(overlap_vars) == 0L) {
      return(data.table::data.table(scenario = scn, status = 'not_run', message = 'Insufficient overlap variables', min_fitted = NA_real_, q25_fitted = NA_real_, median_fitted = NA_real_, mean_fitted = NA_real_, q75_fitted = NA_real_, max_fitted = NA_real_, prop_gt_0_9 = NA_real_, prop_lt_0_1 = NA_real_))
    }
    obs_flag <- data.table::copy(obs[, ..overlap_vars])
    obs_flag[, source := 'observed']
    int_flag <- data.table::copy(scn_dt[, ..overlap_vars])
    int_flag[, source := 'intervention']
    check_dt <- data.table::rbindlist(list(obs_flag, int_flag), fill = TRUE)
    overlap_fit <- tryCatch({
      stats::glm(
        stats::as.formula(paste0("I(source == 'intervention') ~ ", paste(overlap_vars, collapse = ' + '))),
        data = as.data.frame(check_dt),
        family = stats::binomial()
      )
    }, error = function(e) e)
    if (inherits(overlap_fit, 'error')) {
      return(data.table::data.table(scenario = scn, status = 'failed', message = conditionMessage(overlap_fit), min_fitted = NA_real_, q25_fitted = NA_real_, median_fitted = NA_real_, mean_fitted = NA_real_, q75_fitted = NA_real_, max_fitted = NA_real_, prop_gt_0_9 = NA_real_, prop_lt_0_1 = NA_real_))
    }
    fitted_vals <- stats::fitted(overlap_fit)
    q <- stats::quantile(fitted_vals, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE)
    data.table::data.table(
      scenario = scn,
      status = 'ok',
      message = NA_character_,
      min_fitted = q[1],
      q25_fitted = q[2],
      median_fitted = q[3],
      mean_fitted = mean(fitted_vals, na.rm = TRUE),
      q75_fitted = q[4],
      max_fitted = q[5],
      prop_gt_0_9 = mean(fitted_vals > 0.9, na.rm = TRUE),
      prop_lt_0_1 = mean(fitted_vals < 0.1, na.rm = TRUE)
    )
  }), fill = TRUE)

  link_summary <- NULL
  if (!is.null(model_fit) && !is.null(model_fit$outcome_model)) {
    outcome_model <- model_fit$outcome_model
    # speedglm is in Imports so its namespace exists, but its S3 methods
    # (predict.speedglm) are only registered once the namespace is *fully*
    # initialised.  In a fresh R session where speedglm was never actively
    # called, the registration may not have happened yet.
    # loadNamespace() is idempotent: a no-op if already loaded.
    tryCatch(loadNamespace("speedglm"), error = function(e) NULL)
    link_summary <- data.table::rbindlist(lapply(scenarios, function(scn) {
      scn_dt <- int[get(scenario_col) == scn]
      summarize_link <- function(newdata, group_name) {
        pred <- tryCatch({
          p   <- as.numeric(predict_tvcqgcomp_model(outcome_model, newdata = newdata, type = 'response'))
          fam <- tryCatch(family(outcome_model), error = function(e) NULL)
          if (!is.null(fam) && is.function(fam$linkfun)) fam$linkfun(p) else log(p / (1 - p))
        }, error = function(e) e)
        if (inherits(pred, 'error')) {
          return(data.table::data.table(
            scenario = scn,
            group = group_name,
            status = 'failed',
            message = conditionMessage(pred),
            min_link = NA_real_,
            q25_link = NA_real_,
            median_link = NA_real_,
            mean_link = NA_real_,
            q75_link = NA_real_,
            max_link = NA_real_
          ))
        }
        data.table::data.table(
          scenario = scn,
          group = group_name,
          status = 'ok',
          message = NA_character_,
          min_link = min(pred, na.rm = TRUE),
          q25_link = stats::quantile(pred, 0.25, na.rm = TRUE),
          median_link = stats::quantile(pred, 0.5, na.rm = TRUE),
          mean_link = mean(pred, na.rm = TRUE),
          q75_link = stats::quantile(pred, 0.75, na.rm = TRUE),
          max_link = max(pred, na.rm = TRUE)
        )
      }
      data.table::rbindlist(list(
        summarize_link(nat, 'natural'),
        summarize_link(scn_dt, 'intervention')
      ), fill = TRUE)
    }), fill = TRUE)
  }

  list(
    overall_support = overall_support,
    time_support = time_support,
    mean_by_time = mean_by_time,
    joint_support_overall = joint_support_overall,
    joint_support_time = joint_support_time,
    overlap_summary = overlap_summary,
    link_summary = link_summary
  )
}































