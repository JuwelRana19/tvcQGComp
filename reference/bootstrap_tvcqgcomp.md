# Bootstrap tvcQGcomp Fits

Run a subject-level cluster bootstrap for `tvcQGcomp`. In regular mode
the function returns replicate-level meta-model coefficients. In audit
mode it also stores the per-replicate bootstrap data, Monte Carlo data,
natural course output, intervention outputs, and risk curves. The
canonical user-facing survival bootstrap entry point is
`tvcQGComp_survival_boot()`.

## Usage

``` r
tvcQGComp_survival_boot(
  data,
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
  keep_boot_curves = FALSE
)

summarize_tvcqgcomp_bootstrap(boot_res)
```

## Arguments

- data, config, mc_size, intervention_scenarios, natural_course, seed,
  replace_mc, meta_formula, meta_family, meta_target, fit_meta,
  fit_exposure_models, exposure_scale, quantization_breaks,
  joint_effect_fn:

  See
  [`tvcQGComp_survival()`](https://JuwelRana19.github.io/tvcQGComp/reference/run_tvcqgcomp.md).
  When multiple meta targets are requested, bootstrap result columns are
  prefixed by target, for example `coef_HR_joint_effect`,
  `coef_HD_joint_effect`, and `coef_RR_joint_effect`. When
  `fit_meta = FALSE`, replicate outputs contain scenario-specific final
  risks instead of meta-model coefficients.

- n_boot:

  Number of bootstrap replicates.

- parallel:

  Logical; whether to evaluate bootstrap replicates in parallel.

- n_workers:

  Number of workers when `parallel = TRUE`.

- batch_size:

  Number of replicates per parallel batch.

- checkpoint_file:

  Optional RDS file to update during the bootstrap run.

- resume_from_checkpoint:

  Logical; if `TRUE` and `checkpoint_file` exists, resume from the last
  completed bootstrap replicate stored in the checkpoint.

- verbose:

  Logical; print bootstrap progress.

- stop_on_error:

  Logical; if `TRUE`, stop the bootstrap at the first failed replicate
  after saving a checkpoint. If `FALSE`, failed replicates are retained
  with an `error: ...` status and the bootstrap continues.

- audit_mode:

  Logical; if `TRUE`, keep replicate objects for auditing.

- keep_boot_objects:

  Logical; explicit switch to retain audit objects.

- keep_boot_curves:

  Logical; retain a lightweight per-replicate curve table containing
  time-specific risk and survival without storing the full bootstrap
  audit objects.

- boot_res:

  Result from `tvcQGComp_survival_boot()`.

## Value

`tvcQGComp_survival_boot()` returns either a data table of replicate
meta-model coefficients or scenario-specific final risks or, when curve
summaries and/or audit objects are retained, a list with `results` and
optional `curve_results`, `audit_objects`, and `audit_mode`.
`summarize_tvcqgcomp_bootstrap()` returns a summary table with means,
standard errors, and percentile intervals.
