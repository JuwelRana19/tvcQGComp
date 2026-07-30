# Create a tvcQGcomp Configuration and Quantized Scenarios

Create the configuration object used by `tvcQGcomp`, generate simple
static Q1 to QK intervention scenarios, and validate observed data
against the package requirements.

## Usage

``` r
make_tvcqgcomp_config(
  id,
  time_in = NULL,
  time_name = NULL,
  time_out,
  outcome,
  outcome_type = "survival",
  exposures,
  exposure_lags = NULL,
  exposure_formulas = NULL,
  exposure_types = NULL,
  time_varying_covariates,
  tvc_lags = NULL,
  covtypes = NULL,
  passive_tvc_vars = NULL,
  time_fixed_covariates,
  factor_vars = NULL,
  outcome_formula,
  tvc_formulas,
  q = NULL,
  q_levels = NULL,
  natural_course = "calibration",
  auto_history = TRUE,
  baselags = FALSE,
  meta_target = "HR",
  categorical_sim = "stochastic",
  exposome = NULL,
  exposure_scale = NULL,
  bounded_supports = NULL
)

make_q_scenarios(config, levels = NULL, prefix = "Q")

generate_intervention_scenarios(config, levels = NULL, prefix = "Q")

validate_tvcqgcomp_data(data, config)
```

## Arguments

- id:

  Subject identifier column name.

- time_in:

  Start-of-interval time column name.

- time_name:

  Alternative name for `time_in`.

- time_out:

  End-of-interval time column name.

- outcome:

  Outcome column name.

- outcome_type:

  Outcome type. The current public release supports `"survival"` for
  pooled person-time survival analyses.

- exposures:

  Exposure column names.

- exposure_lags:

  Lagged exposure column names. Defaults to
  `paste0(exposures, "_lag1")`.

- exposure_formulas:

  Optional formulas for modeled natural-course exposures.

- exposure_types:

  Named vector of exposure model types. Supported values include
  `"categorical"`, `"multinomial"`, `"binary"`, `"normal"`,
  `"bounded_normal"`, `"bounded_normal_snap"`, and ordinal-logit aliases
  `"ordinal"`, `"ordinal_logit"`, `"ordinal_polr"`, `"ordered_logit"`,
  `"polr"`.

- time_varying_covariates:

  Time-varying confounder column names.

- tvc_lags:

  Lagged time-varying confounder column names.

- covtypes:

  Named vector of time-varying confounder model types. Supported values
  include `"categorical"`, `"multinomial"`, `"binary"`, `"normal"`,
  `"bounded_normal"`, `"bounded_normal_snap"`, and ordinal-logit aliases
  `"ordinal"`, `"ordinal_logit"`, `"ordinal_polr"`, `"ordered_logit"`,
  `"polr"`.

- passive_tvc_vars:

  Optional subset of `time_varying_covariates` to carry forward from the
  sampled observed id-time path without fitting a stochastic TVC model.
  Passive TVCs remain available to outcome, exposure, and other
  confounder formulas at each time point.

- time_fixed_covariates:

  Baseline covariates that stay fixed within subject.

- factor_vars:

  Variables that should be treated as factors during simulation and
  prediction.

- outcome_formula:

  Outcome model formula.

- tvc_formulas:

  Named list of time-varying confounder model formulas.

- q:

  Primary quantization level argument. If `NULL`, the function falls
  back to `q_levels` (or 4).

- q_levels:

  Number of quantile levels used to define the intervention grid.

- natural_course:

  One or more natural-course modes. The preferred names are `"none"`,
  `"calibration"`, `"observed_exposome"`, and `"modeled_exposome"`.
  Backward-compatible aliases `"observed"` and `"modeled"` are also
  accepted and are normalized to `"calibration"` and
  `"modeled_exposome"`, respectively. `"calibration"` predicts outcomes
  on the observed data path without forward simulation;
  `"observed_exposome"` forward-simulates TVCs while keeping sampled
  observed exposure histories fixed; and `"modeled_exposome"`
  forward-simulates both TVCs and exposures.

- auto_history:

  Logical; whether to generate required lag histories from their base
  variables during data preparation. When `TRUE`, lag columns referenced
  by the configuration or model formulas do not need to be present in
  the input data.

- baselags:

  Logical scalar following the `gfoRmula` convention for early lag
  histories when pre-baseline rows are unavailable. If `FALSE`, early
  `*_lag2`, `*_lag3`, etc. values are set to 0 for numeric variables or
  the reference level for factors. If `TRUE`, early lag values are set
  to the subject's baseline value.

- meta_target:

  Second-stage summary target. Use `"HR"` for the pooled hazard-based
  person-time multiplicative meta model, `"HD"` for the pooled
  hazard-based additive meta model, `"RR"` for the end-of-follow-up
  cumulative-risk ratio meta model, or `"RD"` for the end-of-follow-up
  cumulative-risk difference meta model. Multiple targets can be
  requested, for example `c("HR", "HD")` or `c("RR", "RD")`; the legacy
  shortcut `"both"` expands to `c("HR", "RR")`. Historical aliases
  `"hazard"` and `"cumrisk_final"` are also accepted.

- categorical_sim:

  Categorical simulation rule for time-varying covariates and modeled
  exposures. `"stochastic"` uses multinomial probability draws;
  `"class"` uses deterministic class prediction.

- exposome:

  Exposure representation. Use `"quantized"` for the current public
  release.

- exposure_scale:

  Backward-compatible alias for `exposome`.

- bounded_supports:

  Optional named list specifying legal support grids for bounded-normal
  snap/ordinal variables. Each entry can be a numeric vector (for
  example `1:5`), `list(values = ...)`, `list(min=, max=, by=)`, or
  `"observed"` to use observed unique values.

- data:

  Observed data to validate.

- config:

  A `tvcQGcomp` configuration object.

- levels:

  Intervention levels to assign to all exposures. Defaults to
  `0:(q_levels - 1)`.

- prefix:

  Scenario name prefix.

## Value

`make_tvcqgcomp_config()` returns a named list used by the package.
`make_q_scenarios()` and `generate_intervention_scenarios()` return a
named list of static intervention vectors. `validate_tvcqgcomp_data()`
returns `TRUE` invisibly and throws an error if the input data do not
satisfy the required structure.
