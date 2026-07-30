test_that("parallel bootstrap argument validation rejects bad worker settings", {
  expect_error(
    validate_bootstrap_parallel_args(
      parallel = TRUE,
      n_workers = 0,
      call_name = "bootstrap_tvcqgcomp"
    ),
    "n_workers"
  )

  expect_error(
    validate_bootstrap_parallel_args(
      parallel = TRUE,
      n_workers = NA,
      call_name = "bootstrap_tvcqgcomp"
    ),
    "n_workers"
  )

  expect_error(
    validate_bootstrap_parallel_args(
      parallel = TRUE,
      n_workers = 2,
      batch_size = 0,
      call_name = "bootstrap_tvcqgcomp"
    ),
    "batch_size"
  )
})

test_that("parallel bootstrap argument validation warns on retained objects", {
  expect_warning(
    validate_bootstrap_parallel_args(
      parallel = TRUE,
      n_workers = 2,
      audit_mode = TRUE,
      call_name = "bootstrap_tvcqgcomp"
    ),
    "substantial RAM"
  )
})
