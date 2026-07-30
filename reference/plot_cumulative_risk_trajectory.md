# Plot Cumulative Risk Trajectories Over Time

Draw cumulative-risk trajectories over time for natural-course and
intervention scenarios from a `tvcQGcomp` run object or a pre-summarized
trajectory table. This is the primary plotting function for cumulative
incidence.

## Usage

``` r
plot_cumulative_risk_trajectory(
  x,
  include_natural = TRUE,
  natural_label = "Natural Course",
  include_average = FALSE,
  average_scenarios = paste0("Q", 0:3),
  average_label = "Average (Q0-Q3)",
  average_col = "black",
  scenario_col = "scenario",
  time_col = "TimeInn",
  risk_col = "mean_risk",
  percent = TRUE,
  scenario_order = NULL,
  scenario_labels = NULL,
  cols = NULL,
  include_points = TRUE,
  pch = 16,
  lwd = 2,
  xlab = "Follow-up time",
  ylab = if (percent) "Cumulative Risk (pct)" else "Cumulative Risk",
  main = "Cumulative risk over time",
  legend_pos = "topleft",
  ...
)
```

## Arguments

- x:

  A `tvcQGcomp` run object or a data.frame/data.table containing
  summarized risk trajectories.

- include_natural:

  Logical; when `x` is a run object, include the natural-course curve if
  available.

- natural_label:

  Label to use for the natural-course curve in the legend.

- include_average:

  Logical; add a line for the average of `average_scenarios` at each
  follow-up time.

- average_scenarios:

  Character vector of scenario names to average. Defaults to `Q0`
  through `Q3`.

- average_label:

  Legend label for the averaged trajectory.

- average_col:

  Color for the averaged trajectory when `cols` is not supplied.

- scenario_col:

  Name of the scenario column.

- time_col:

  Name of the follow-up time column.

- risk_col:

  Name of the cumulative risk column to plot.

- percent:

  Logical; if `TRUE`, plot risk on the percent scale.

- scenario_order:

  Optional character vector giving the desired plotting order of
  scenarios.

- scenario_labels:

  Optional replacement labels for scenarios. May be a character vector
  named by scenario or an unnamed vector with one label per scenario.

- cols:

  Optional vector of plotting colors. May be named by scenario.

- include_points:

  Logical; whether to add points at each time value.

- pch:

  Point character used when `include_points = TRUE`.

- lwd:

  Line width.

- xlab:

  X-axis label.

- ylab:

  Y-axis label.

- main:

  Plot title.

- legend_pos:

  Legend position passed to
  [`legend`](https://rdrr.io/r/graphics/legend.html).

- ...:

  Additional arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

The plotting helper returns the plotted data invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- tvcQGComp_survival(NO2_dt, cfg, mc_size = 100, seed = 1234)
plot_cumulative_risk_trajectory(res)
plot_cumulative_risk_trajectory(res, include_natural = FALSE, include_average = TRUE)
} # }
```
