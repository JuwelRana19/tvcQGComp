# Diagnose Drift Between Observed and Natural-Course TVCs

Compare observed and simulated natural-course trajectories over time for
time-varying covariates. Categorical variables are summarized through
time-specific level proportions, while continuous variables are
summarized through time-specific means, standard deviations, and
selected quantiles.

## Usage

``` r
diagnose_tvcqgcomp_drift(observed_data = NULL, natural = NULL, run_obj = NULL,
  config = NULL, variables = NULL, include = c("both", "categorical",
  "continuous"), continuous_probs = c(0.1, 0.5, 0.9))
```

## Arguments

- observed_data:

  Observed longitudinal data used as the calibration reference.

- natural:

  Simulated natural-course data to compare against the observed data.

- run_obj:

  Optional object returned by
  [`tvcQGComp_survival()`](https://JuwelRana19.github.io/tvcQGComp/reference/run_tvcqgcomp.md).
  When supplied, `config`, `observed_data`, and `natural` are filled
  from the run object when available.

- config:

  Optional configuration object. If omitted, `run_obj$config` is used.

- variables:

  Optional character vector of variables to check. Defaults to the
  configured time-varying covariates.

- include:

  Whether to return checks for `"categorical"`, `"continuous"`, or
  `"both"` types.

- continuous_probs:

  Quantiles to compare over time for continuous variables.

## Value

A list with the inferred variable types and separate detailed and
summary tables for categorical and continuous drift over time.
