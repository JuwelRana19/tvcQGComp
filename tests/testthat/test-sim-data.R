test_that("bundled simulation data use continuous exposures without lag columns", {
  data("sim_data", package = "tvcQGComp")
  exposures <- c("BC", "NIT", "SO4", "NH4", "OM")

  expect_equal(nrow(sim_data), 5584L)
  expect_equal(ncol(sim_data), 28L)
  expect_equal(length(unique(sim_data$UniqID)), 500L)
  expect_equal(sum(sim_data$status), 100)
  expect_false(any(grepl("_lag[0-9]+$", names(sim_data))))

  for (exposure in exposures) {
    values <- sim_data[[exposure]]
    expect_true(is.numeric(values))
    expect_true(all(is.finite(values)))
    expect_gt(length(unique(values)), 100L)
    expect_true(any(abs(values - round(values)) > 1e-8))
  }
})

test_that("automatic history creates lag-one variables absent from input data", {
  config <- list(
    id = "id",
    time_in = "time",
    exposures = "x",
    exposure_lags = "x_lag1",
    tvc = "z",
    tvc_lags = "z_lag1",
    auto_history = TRUE,
    baselags = TRUE,
    outcome_formula = y ~ x + x_lag1 + z + z_lag1,
    tvc_formulas = list(z = z ~ x_lag1),
    exposure_formulas = NULL
  )
  input <- data.table::data.table(
    id = rep(1:2, each = 3),
    time = rep(0:2, times = 2),
    x = 1:6,
    z = factor(c("a", "b", "a", "b", "a", "b")),
    y = 0
  )

  result <- generate_history_columns(input, config)

  expect_true(all(c("x_lag1", "z_lag1") %in% names(result)))
  expect_equal(result$x_lag1, c(1L, 1L, 2L, 4L, 4L, 5L))
  expect_equal(
    as.character(result$z_lag1),
    c("a", "a", "b", "b", "b", "a")
  )
})
