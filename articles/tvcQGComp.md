# Get started with tvcQGComp

## Overview

`tvcQGComp` estimates the survival effect of joint interventions on a
quantized time-varying exposure mixture. The core workflow is:

1.  arrange the observed data in long person-period format;
2.  define the outcome, exposures, covariates, histories, and estimand;
3.  validate the data against that configuration;
4.  generate quantized intervention scenarios;
5.  run Monte Carlo g-computation; and
6.  use a subject-level bootstrap for confidence intervals.

This public release focuses on main-effect survival analyses with
quantized exposures.

## Required data structure

The input must contain one row per subject and follow-up interval. At
minimum, the configured columns must identify:

| Role | Example |
|----|----|
| Subject identifier | `UniqID` |
| Interval start and end | `TimeInn`, `TimeOut` |
| Interval event indicator | `status` |
| Time-varying exposures | `BC`, `NIT`, `SO4`, `NH4`, `OM` |
| Time-varying covariates | `CSize`, `Urban_form`, `CanadianRegion`, neighborhood indices |
| Baseline covariates | `sex`, education, employment, and other demographic variables |

Subject identifiers and interval times should uniquely identify rows.
Base variables used by model formulas must be present, correctly typed,
and supported at every time required by the simulation. Lagged formula
terms do not need to be stored in the input when `auto_history = TRUE`.

## Configure the analysis

The package includes the complete simulated dataset used in this
example.

``` r

library(tvcQGComp)

data("sim_data", package = "tvcQGComp")
dat <- sim_data
dat$TimeOut <- as.integer(as.character(dat$TimeOut))

exposures <- c("BC", "NIT", "SO4", "NH4", "OM")
active_tvc <- c(
  "Urban_form", "CanadianRegion", "CSize", "deprivation",
  "dependency", "instability", "ethnicconcentration"
)
baseline <- c(
  "age", "income_inadequacy", "sex", "visible_minority",
  "IndigenousIdentity", "marsth", "Education", "employment", "Occupation"
)

# These are formula term names only. auto_history creates the columns.
active_lag_terms <- paste0(active_tvc, "_lag1")
tvc_formulas <- setNames(lapply(active_tvc, function(response) {
  reformulate(
    c(
      "TimeOut", "I(TimeOut^2)", baseline,
      setdiff(active_lag_terms, paste0(response, "_lag1"))
    ),
    response = response
  )
}), active_tvc)

cfg <- make_tvcqgcomp_config(
  id = "UniqID",
  time_in = "TimeInn",
  time_out = "TimeOut",
  outcome = "status",
  exposures = exposures,
  time_varying_covariates = active_tvc,
  time_fixed_covariates = baseline,
  factor_vars = c("income_inadequacy", active_tvc, baseline),
  outcome_formula = reformulate(
    c(exposures, "TimeOut", "I(TimeOut^2)", baseline, active_tvc),
    response = "status"
  ),
  tvc_formulas = tvc_formulas,
  covtypes = c(
    setNames(rep("categorical", length(active_tvc)), active_tvc)
  ),
  q = 4,
  exposome = "quantized",
  natural_course = "calibration",
  auto_history = TRUE,
  baselags = TRUE,
  meta_target = c("HR", "RR", "RD")
)

validate_tvcqgcomp_data(dat, cfg)
```

The configuration is the analysis contract. It should be finalized
before fitting the point estimate or bootstrap so both stages use
identical models and interventions. In this simple example, `age` and
`income_inadequacy` are used as baseline/time-fixed adjustment
variables, while the seven neighborhood and urban-form covariates are
simulated forward as active time-varying covariates.

## Define interventions

The input exposures are continuous. The package determines their
quartile cutpoints internally, and the default helper assigns every
internally quantized exposure to each level in turn.

``` r

scenarios <- generate_intervention_scenarios(cfg)
names(scenarios)
```

Inspect `scenarios` before fitting. The scenario names become labels in
the result tables and plots.

## Fit the point estimate

``` r

fit <- tvcQGComp_survival(
  data = dat,
  config = cfg,
  mc_size = 10000,
  intervention_scenarios = scenarios,
  seed = 1234,
  verbose = TRUE
)

fit$risk_trajectory
fit$meta_effect_summary
```

The returned object contains fitted nuisance models, simulated
intervention data, risk summaries, natural-course outputs when
requested, and the pooled second-stage model.

## Inspect risk trajectories

``` r

risk_curves <- as_risk_curve_df(fit)
head(risk_curves)

plot_cumulative_risk_trajectory(fit)
```

Risk trajectories are useful for understanding when scenarios begin to
separate. Final HR, RR, and RD estimates come from their configured
second-stage meta-models.

## Next steps

Use the bootstrap for confidence intervals and inspect diagnostics
before reporting results. See the [bootstrap and
diagnostics](https://JuwelRana19.github.io/tvcQGComp/articles/bootstrap-and-diagnostics.md)
guide.
