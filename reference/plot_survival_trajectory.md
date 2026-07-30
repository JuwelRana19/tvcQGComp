# Plot Survival Trajectories Over Time

Draw survival trajectories over time for natural-course and intervention
scenarios from a `tvcQGcomp` run object or a pre-summarized trajectory
table.

## Usage

``` r
plot_survival_trajectory(
  x,
  include_natural = TRUE,
  natural_label = "Natural Course",
  include_average = FALSE,
  average_scenarios = paste0("Q", 0:3),
  average_label = "Average (Q0-Q3)",
  average_col = "black",
  scenario_col = "scenario",
  time_col = "TimeInn",
  survival_col = "mean_survival",
  percent = TRUE,
  scenario_order = NULL,
  scenario_labels = NULL,
  cols = NULL,
  include_points = TRUE,
  pch = 16,
  lwd = 2,
  xlab = "Follow-up time",
  ylab = if (percent) "Survival (%)" else "Survival probability",
  main = "Survival over time",
  legend_pos = "bottomleft",
  ...
)
```

## Arguments

- x:

  A `tvcQGcomp` run object or a data frame containing summarized
  trajectories.

- include_natural:

  Logical; include the natural-course curve when available.

- natural_label:

  Legend label for the natural-course curve.

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

- survival_col:

  Name of the survival column to plot.

- percent:

  Logical; plot survival on the percentage scale.

- scenario_order:

  Optional character vector giving scenario order.

- scenario_labels:

  Optional replacement labels for scenarios.

- cols:

  Optional plotting colors, optionally named by scenario.

- include_points:

  Logical; add points at each follow-up time.

- pch:

  Point character.

- lwd:

  Line width.

- xlab:

  Horizontal-axis label.

- ylab:

  Vertical-axis label.

- main:

  Plot title.

- legend_pos:

  Legend position passed to
  [`legend`](https://rdrr.io/r/graphics/legend.html).

- ...:

  Additional arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Returns the plotted trajectory data invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- tvcQGComp_survival(dat, cfg, intervention_scenarios = scenarios)
plot_survival_trajectory(fit, include_natural = FALSE)
plot_survival_trajectory(fit, include_natural = FALSE, include_average = TRUE)
} # }
```
