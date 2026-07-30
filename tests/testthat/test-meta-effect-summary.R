test_that("direct contrast configuration and helpers are absent", {
  obsolete_args <- c("direct_contrasts", "final_risk_reference")
  obsolete_helpers <- c(
    "resolve_curve_reference_scenario",
    "add_curve_contrasts"
  )
  namespace <- asNamespace("tvcQGComp")

  expect_false(any(obsolete_args %in% names(formals(make_tvcqgcomp_config))))
  expect_false(any(vapply(
    obsolete_helpers,
    exists,
    logical(1),
    envir = namespace,
    inherits = FALSE
  )))
})

test_that("bootstrap endpoint extraction is derived from risk trajectories", {
  trajectory <- data.table::data.table(
    scenario = rep(c("Q0", "Q1"), each = 2),
    TimeInn = rep(0:1, times = 2),
    TimeOut = rep(1:2, times = 2),
    mean_risk = c(0.05, 0.10, 0.08, 0.15),
    mean_survival = c(0.95, 0.90, 0.92, 0.85)
  )
  config <- list(time_in = "TimeInn", time_out = "TimeOut")

  result <- flatten_risk_trajectory_endpoint(trajectory, config)
  expect_named(result, c("risk_Q0", "risk_Q1"))
  expect_equal(unname(unlist(result)), c(0.10, 0.15))
})

test_that("RR and RD summaries use their meta-model joint-effect coefficients", {
  meta_data <- data.frame(
    cumrisk = c(0.10, 0.11, 0.13, 0.14, 0.17, 0.18, 0.21, 0.22),
    joint_effect = rep(0:3, each = 2)
  )
  rr_model <- stats::glm(
    cumrisk ~ joint_effect,
    data = meta_data,
    family = stats::quasibinomial("log")
  )
  rd_model <- stats::glm(
    cumrisk ~ joint_effect,
    data = meta_data,
    family = stats::quasibinomial("identity")
  )

  result <- summarize_meta_model_effects(list(RR = rr_model, RD = rd_model))

  expect_equal(
    result[meta_target == "RR", estimate],
    exp(stats::coef(rr_model)[["joint_effect"]])
  )
  expect_equal(
    result[meta_target == "RD", estimate],
    unname(stats::coef(rd_model)[["joint_effect"]])
  )
  expect_named(result, c("meta_target", "estimate", "increment"))
  expect_equal(result$increment, c(1, 1))
})

test_that("run printing shows a concise meta-model effect summary", {
  meta_effect_summary <- data.table::data.table(
    meta_target = c("HR", "RR", "RD"),
    estimate = c(exp(0.20), exp(0.18), 0.04),
    increment = 1
  )
  run <- structure(
    list(
      config = list(
        outcome_type = "survival",
        meta_target = c("HR", "RR", "RD"),
        categorical_sim = "stochastic",
        natural_course = "none"
      ),
      interventions = list(Q0 = 0, Q1 = 1),
      mc_data = data.frame(id = 1:2),
      risk_trajectory = data.table::data.table(
        scenario = rep(c("Q0", "Q1"), each = 2),
        TimeInn = rep(0:1, times = 2),
        mean_risk = c(0, 0.10, 0, 0.15),
        mean_survival = c(1, 0.90, 1, 0.85)
      ),
      meta_effect_summary = meta_effect_summary
    ),
    class = "tvcqgcomp_run"
  )

  printed <- capture.output(print(run))

  expect_false(any(grepl("prediction summary", printed, fixed = TRUE)))
  expect_false(any(grepl("Mean Survival", printed, fixed = TRUE)))
  trajectory_printed <- capture.output(print(run, prediction = "trajectory"))
  expect_true(any(grepl("Risk trajectory", trajectory_printed, fixed = TRUE)))
  expect_true(any(grepl("Scenario", trajectory_printed, fixed = TRUE)))
  expect_true(any(grepl("Mean Risk", trajectory_printed, fixed = TRUE)))
  expect_true(any(grepl("Mean Survival", trajectory_printed, fixed = TRUE)))
  expect_true(any(grepl("Meta Target", printed, fixed = TRUE)))
  expect_true(any(grepl("Estimate", printed, fixed = TRUE)))
  expect_true(any(grepl("Increment", printed, fixed = TRUE)))
  expect_false(any(grepl("effect_measure|model_link|link_estimate|joint_effect", printed)))
})

test_that("public plotting uses cumulative risk with optional survival", {
  run <- structure(
    list(
      config = list(time_in = "TimeInn"),
      interventions = list(Q0 = 0, Q1 = 1, Q2 = 2, Q3 = 3),
      risk_trajectory = data.table::data.table(
        scenario = rep(c("Q0", "Q1", "Q2", "Q3"), each = 2),
        TimeInn = rep(0:1, times = 4),
        mean_risk = c(0, 0.10, 0, 0.15, 0, 0.20, 0, 0.25),
        mean_survival = c(1, 0.90, 1, 0.85, 1, 0.80, 1, 0.75)
      )
    ),
    class = c("tvcqgcomp_run", "list")
  )
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  risk_data <- plot(run, include_natural = FALSE, include_average = TRUE)
  survival_data <- plot_survival_trajectory(run, include_natural = FALSE, include_average = TRUE)

  expect_s3_class(risk_data, "data.table")
  expect_s3_class(survival_data, "data.table")
  avg_risk <- risk_data[scenario == "Average (Q0-Q3)"]
  avg_survival <- survival_data[scenario == "Average (Q0-Q3)"]
  expect_equal(avg_risk$mean_risk, c(0, 0.175))
  expect_equal(avg_survival$mean_survival, c(1, 0.825))
})
