# Run the tvcQGcomp Pipeline

Fit the nuisance models, create the Monte Carlo skeleton, simulate the
natural course and static quantized interventions, and fit the pooled
person-time meta model. The canonical user-facing survival entry point
is `tvcQGComp_survival()`.

## Usage

``` r
tvcQGComp_survival(
  data,
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
  verbose = FALSE
)
```

## Arguments

- data:

  Observed input data.

- config:

  A `tvcQGcomp` configuration object.

- mc_size:

  Monte Carlo sample size.

- intervention_scenarios:

  Named list of static quantized intervention values.

- natural_course:

  Any subset of the preferred labels `"calibration"`,
  `"observed_exposome"`, and `"modeled_exposome"`, or `"none"`.
  Backward-compatible aliases `"observed"` and `"modeled"` are also
  accepted. `"calibration"` predicts outcomes on the observed data path,
  `"observed_exposome"` forward-simulates TVCs while keeping sampled
  observed exposure histories, and `"modeled_exposome"`
  forward-simulates both TVCs and exposures.

- seed:

  Random seed for Monte Carlo resampling and simulation.

- replace_mc:

  Whether to sample subjects with replacement when creating the Monte
  Carlo skeleton.

- meta_formula:

  Optional pooled meta-model formula.

- meta_family:

  Optional model family for the pooled meta model. If a single target is
  requested and `meta_family` is `NULL`, defaults are
  `quasibinomial("logit")` for `"HR"`, `gaussian("identity")` for
  `"HD"`, `quasibinomial("identity")` for `"RD"`, and
  `quasibinomial("log")` for `"RR"`. When multiple meta targets are
  requested, pass either `NULL` to use defaults for each target or a
  named list such as
  `list(HR = quasibinomial("logit"), HD = gaussian("identity"))`.

- meta_target:

  Second-stage summary target. Use `"HR"`, `"HD"`, `"RR"`, `"RD"`, or a
  vector such as `c("HR", "HD")`. Historical aliases `"hazard"` and
  `"cumrisk_final"` are also accepted; the legacy shortcut `"both"`
  expands to `c("HR", "RR")`. When multiple targets are requested, the
  first-stage simulation is run once and all requested second-stage
  meta-models are fit from the same simulated intervention data.

- fit_meta:

  Optional logical controlling whether the second-stage meta-model is
  fit. Defaults to `TRUE` for quantized exposures. When `FALSE`,
  `risk_trajectory` still reports time-specific scenario risks and
  survival, but no HR, HD, RR, or RD effect estimate is produced.

- fit_exposure_models:

  Optional logical controlling whether exposure evolution models are
  fit. Defaults to `FALSE` when `natural_course = "none"` and every
  intervention scenario deterministically assigns every exposure;
  otherwise defaults to `TRUE`.

- exposure_scale:

  Exposure representation used during first-stage fitting and
  simulation. For this minimal release, use `"quantized"`.

- quantization_breaks:

  Optional fixed cutpoints used when `exposure_scale = "quantized"`.

- joint_effect_fn:

  Optional function mapping a scenario vector to a scalar joint effect.

- verbose:

  Logical; print progress and status messages during fitting and
  simulation.

## Value

A list containing `mc_data`, `natural`, `interventions`,
`intervention_data`, `risk_trajectory`, `meta_model`, `meta_models`,
`meta_effect_summary`, calibration outputs, and the fitted nuisance
models. `risk_trajectory` is the sole public scenario-prediction table
and reports time-specific cumulative risk and survival by intervention
scenario. HR, HD, RR, and RD estimates are derived exclusively from the
requested second-stage meta-models and reported in
`meta_effect_summary`, which contains only `meta_target`, `estimate`,
and `increment`. For HR and RR, `estimate` is the exponentiated
`joint_effect` coefficient; for HD and RD, it is the untransformed
coefficient. Each estimate represents a one-unit increase in the joint
quantile intervention. When `fit_meta = FALSE`, the meta-model elements
and `meta_effect_summary` are `NULL`. When `natural_course = "none"`,
`natural` and the natural-course calibration outputs are `NULL`. The
Monte Carlo skeleton inherits baseline and lagged predictors from the
observed data, so simulated baseline/history support must be complete
for the predictors required by the nuisance models. The `calibration`
element includes `calibration`/`calibration_curve`,
`observed_exposome`/`observed_exposome_curve`, and
`modeled_exposome`/`modeled_exposome_curve`; backward-compatible aliases
`observed` and `modeled` are retained.
