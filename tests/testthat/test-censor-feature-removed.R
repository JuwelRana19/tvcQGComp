test_that("censor-model configuration arguments are not public", {
  censor_args <- c(
    "censor_name",
    "censor_formula",
    "censor_mode",
    "ipw_cutoff_quantile",
    "ipw_cutoff_value"
  )

  expect_false(any(censor_args %in% names(formals(make_tvcqgcomp_config))))
})

test_that("censor-model helpers are absent from the internal engine", {
  censor_helpers <- c(
    "has_explicit_censoring",
    "use_simulated_censoring",
    "compute_ipcw_weights",
    "compute_weighted_observed_survival"
  )
  namespace <- asNamespace("tvcQGComp")

  expect_false(any(vapply(
    censor_helpers,
    exists,
    logical(1),
    envir = namespace,
    inherits = FALSE
  )))
})
