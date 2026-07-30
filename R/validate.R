get_formula_var_names <- function(config) {
  tvc_formulas <- Filter(Negate(is.null), config$tvc_formulas)
  exposure_formulas <- Filter(Negate(is.null), config$exposure_formulas %||% list())
  unique(c(
    all.vars(config$outcome_formula),
    if (length(tvc_formulas)) unlist(lapply(tvc_formulas, all.vars), use.names = FALSE) else character(),
    if (length(exposure_formulas)) unlist(lapply(exposure_formulas, all.vars), use.names = FALSE) else character()
  ))
}

extract_lag_specs <- function(vars) {
  vars <- unique(as.character(vars))
  out <- lapply(vars, function(v) {
    m <- regexec("^(.*)_lag([0-9]+)$", v)
    reg <- regmatches(v, m)[[1L]]
    if (length(reg) != 3L) {
      return(NULL)
    }
    list(
      lag_name = v,
      base_name = reg[[2L]],
      lag_n = as.integer(reg[[3L]])
    )
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) {
    return(data.table::data.table(
      lag_name = character(),
      base_name = character(),
      lag_n = integer()
    ))
  }
  specs <- data.table::rbindlist(out)
  unique(specs[order(lag_n, lag_name)])
}

get_history_plan <- function(config) {
  formula_vars <- get_formula_var_names(config)
  configured_lags <- c(config$exposure_lags, config$tvc_lags)
  lag_specs <- extract_lag_specs(c(configured_lags, formula_vars))
  default_bases <- stats::setNames(c(config$exposures, config$tvc), c(config$exposure_lags, config$tvc_lags))
  if (nrow(lag_specs)) {
    lag_specs[lag_name %in% names(default_bases), base_name := default_bases[lag_name]]
  }
  list(
    lag_map = stats::setNames(lag_specs$base_name, lag_specs$lag_name),
    lag_specs = lag_specs,
    formula_vars = formula_vars
  )
}

generate_history_columns <- function(dt, config) {
  dt <- data.table::as.data.table(data.table::copy(dt))
  data.table::setorderv(dt, c(config$id, config$time_in))

  if (!isTRUE(config$auto_history)) {
    return(dt)
  }

  plan <- get_history_plan(config)
  lag_specs <- plan$lag_specs
  for (i in seq_len(nrow(lag_specs))) {
    lag_name <- lag_specs$lag_name[[i]]
    base_name <- lag_specs$base_name[[i]]
    lag_n <- lag_specs$lag_n[[i]]
    if (!(base_name %in% names(dt))) next

    if (!(lag_name %in% names(dt))) {
      dt[, (lag_name) := data.table::shift(get(base_name), n = lag_n, type = "lag"), by = c(config$id)]
    }

    if (is.factor(dt[[base_name]])) {
      dt[, (lag_name) := factor(as.character(get(lag_name)), levels = levels(dt[[base_name]]), ordered = is.ordered(dt[[base_name]]))]
    } else if (is.integer(dt[[base_name]])) {
      dt[, (lag_name) := as.integer(get(lag_name))]
    } else if (is.numeric(dt[[base_name]])) {
      dt[, (lag_name) := as.numeric(get(lag_name))]
    }

    if (isTRUE(config$baselags)) {
      dt[, (lag_name) := {
        lag_values <- .SD[[lag_name]]
        reference_value <- .SD[[base_name]][[1L]]
        lag_values[is.na(lag_values)] <- reference_value
        lag_values
      }, by = c(config$id), .SDcols = c(lag_name, base_name)]
    }
  }

  dt
}

prepare_tvcqgcomp_data <- function(data, config) {
  dt <- data.table::as.data.table(data.table::copy(data))

  if (anyDuplicated(names(dt))) {
    stop("Data has duplicated column names.")
  }

  if (!all(c(config$id, config$time_in) %in% names(dt))) {
    stop("Data must include id and time columns before preparation.")
  }

  dt <- generate_history_columns(dt, config)
  data.table::setorderv(dt, c(config$id, config$time_in))

  required <- unique(c(
    config$id,
    config$time_in,
    config$time_out,
    config$outcome,
    config$exposures,
    config$exposure_lags,
    config$tvc,
    config$tvc_lags,
    config$baseline_covariates,
    get_formula_var_names(config)
  ))

  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  dup_n <- dt[, .N, by = c(config$id, config$time_in)][N > 1L, .N]
  if (dup_n > 0L) {
    stop("Data has duplicate id-time rows. Resolve duplicates before fitting or simulation.")
  }

  history_plan <- get_history_plan(config)
  history_cols <- if (nrow(history_plan$lag_specs)) history_plan$lag_specs$lag_name else character()
  non_history_required <- setdiff(required, history_cols)
  check_cols <- intersect(non_history_required, names(dt))
  missing_non_history <- dt[, rowSums(is.na(.SD)) > 0L, .SDcols = check_cols]
  if (any(missing_non_history)) {
    stop("Data contains missing values in required non-history variables. Please use complete-case data for tvcQGcomp 0.2.0.")
  }

  dt
}

validate_tvcqgcomp_data <- function(data, config) {
  prepare_tvcqgcomp_data(data, config)
  invisible(TRUE)
}


