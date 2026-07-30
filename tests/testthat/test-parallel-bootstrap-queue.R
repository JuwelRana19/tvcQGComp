test_that("parallel bootstrap chunking keeps at least n_workers jobs available", {
  chunks <- split_parallel_bootstrap_ids(
    remaining_ids = 1:10,
    n_workers = 4,
    batch_size = 1
  )

  expect_length(chunks[[1L]], 4)
  expect_equal(unlist(chunks, use.names = FALSE), 1:10)
})

test_that("parallel bootstrap chunking uses whole queue when batch_size is NULL", {
  chunks <- split_parallel_bootstrap_ids(
    remaining_ids = 1:6,
    n_workers = 3,
    batch_size = NULL
  )

  expect_length(chunks, 1L)
  expect_length(chunks[[1L]], 6L)
})
