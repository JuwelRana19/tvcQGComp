extract_factor_wrapped_vars <- function(formula) {
  txt <- paste(deparse(formula), collapse = "")
  hits <- gregexpr("(?:as\\.)?factor\\(([^,\\)]+)", txt, perl = TRUE)
  vals <- regmatches(txt, hits)[[1L]]
  if (!length(vals) || identical(vals[[1L]], "-1")) {
    return(character())
  }
  unique(trimws(sub("^(?:as\\.)?factor\\(([^,\\)]+).*$", "\\1", vals, perl = TRUE)))
}

capture_factor_prediction_support <- function(formula, data) {
  vars <- extract_factor_wrapped_vars(formula)
  vars <- intersect(vars, names(data))
  if (!length(vars)) {
    return(NULL)
  }

  out <- lapply(vars, function(v) {
    x <- data[[v]]
    if (is.factor(x)) {
      list(var = v, kind = "factor", levels = levels(x), ordered = is.ordered(x))
    } else if (is.character(x)) {
      list(var = v, kind = "character", levels = sort(unique(stats::na.omit(x))))
    } else if (is.integer(x)) {
      list(var = v, kind = "integer", support = sort(unique(stats::na.omit(x))))
    } else if (is.numeric(x)) {
      list(var = v, kind = "numeric", support = sort(unique(stats::na.omit(x))))
    } else {
      NULL
    }
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) {
    return(NULL)
  }
  names(out) <- vapply(out, `[[`, character(1), "var")
  out
}

attach_prediction_support <- function(model, formula, data) {
  support <- capture_factor_prediction_support(formula, data)
  if (is.null(support)) {
    return(model)
  }
  attr(model, "tvcqgcomp_factor_support") <- support
  if (inherits(model, "tvcqgcomp_bounded_normal")) {
    attr(model$model, "tvcqgcomp_factor_support") <- support
  }
  model
}

nearest_observed_support <- function(x, support) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & !is.na(x)
  support <- sort(unique(stats::na.omit(as.numeric(support))))
  if (!length(support)) {
    return(out)
  }
  if (length(support) == 1L) {
    out[ok] <- support[[1L]]
    return(out)
  }

  xv <- as.numeric(x[ok])
  idx <- findInterval(xv, support)
  idx[idx < 1L] <- 1L
  idx[idx >= length(support)] <- length(support) - 1L
  left <- support[idx]
  right <- support[idx + 1L]
  use_right <- abs(xv - right) < abs(xv - left)
  snapped <- left
  snapped[use_right] <- right[use_right]
  out[ok] <- snapped
  out
}

prepare_newdata_for_prediction <- function(model, newdata) {
  support <- attr(model, "tvcqgcomp_factor_support")
  if (is.null(support) || !length(support)) {
    return(newdata)
  }

  out <- data.table::copy(newdata)
  for (spec in support) {
    v <- spec$var
    if (!(v %in% names(out))) next
    if (identical(spec$kind, "factor")) {
      out[[v]] <- factor(as.character(out[[v]]), levels = spec$levels, ordered = isTRUE(spec$ordered))
    } else if (identical(spec$kind, "character")) {
      out[[v]] <- as.character(out[[v]])
    } else if (identical(spec$kind, "integer")) {
      out[[v]] <- as.integer(round(nearest_observed_support(out[[v]], spec$support)))
    } else if (identical(spec$kind, "numeric")) {
      out[[v]] <- nearest_observed_support(out[[v]], spec$support)
    }
  }
  out
}

predict_tvcqgcomp_model <- function(model, newdata, type = "response", ...) {
  pred_data <- prepare_newdata_for_prediction(model, newdata)
  if (inherits(model, "tvcqgcomp_bounded_normal")) {
    return(stats::predict(model$model, newdata = pred_data, type = type, ...))
  }
  if (inherits(model, c("speedglm", "speedlm"))) {
    # Ensure S3 predict methods are registered in fresh sessions before dispatch.
    tryCatch(loadNamespace("speedglm"), error = function(e) NULL)
    pred <- stats::predict(model, newdata = pred_data, type = type, ...)
    if (identical(type, "response")) {
      fam <- tryCatch(stats::family(model), error = function(e) NULL)
      if (!is.null(fam) && is.function(fam$linkinv)) {
        finite_pred <- suppressWarnings(as.numeric(pred))
        if (length(finite_pred) &&
            any(is.finite(finite_pred)) &&
            any(finite_pred < 0 | finite_pred > 1, na.rm = TRUE)) {
          pred <- fam$linkinv(finite_pred)
        }
      }
    }
    return(pred)
  }
  stats::predict(model, newdata = pred_data, type = type, ...)
}
