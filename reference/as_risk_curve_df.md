# Extract Tidy Risk-Curve Summaries

Combine natural-course and intervention risk-trajectory summaries from a
`tvcQGcomp` run object into one tidy table for downstream plotting or
tabulation.

## Usage

``` r
as_risk_curve_df(
  run_obj,
  include_natural = TRUE,
  natural_label = "Natural Course",
  time_col = NULL,
  scenario_col = "scenario"
)
```

## Arguments

- run_obj:

  A `tvcQGcomp` run object containing `natural_curve` and/or
  `risk_trajectory`. Older objects using `intervention_curve` remain
  supported.

- include_natural:

  Logical; include the natural-course summary if available.

- natural_label:

  Scenario label to assign to the natural-course rows.

- time_col:

  Optional time column name. Defaults to `run_obj$config$time_in`.

- scenario_col:

  Name of the scenario column in the returned table.

## Value

A `data.table` containing one row per scenario and time point, with
columns such as `mean_risk`, `mean_risk_per_1000`, `mean_survival`, and
`mean_survival_per_1000`.

## Details

This helper is useful when you want to reproduce a `gfoRmula`-style
risk-curve plot using `ggplot2` or custom tabulation code.

## Examples

``` r
if (FALSE) { # \dontrun{
curve_dt <- as_risk_curve_df(res)

ggplot(curve_dt, aes(x = TimeInn, y = mean_risk, color = scenario)) +
  geom_line()
} # }
```
