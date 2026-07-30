# Plot Observed vs Natural-Course Drift Over Time

Create small-multiple calibration plots that compare observed
longitudinal trajectories with predicted natural-course trajectories
from `tvcQGcomp`. Continuous variables can be shown using means,
standard deviations, or selected quantiles. Categorical variables can be
shown using either time-specific mean scores or the proportion in a
selected level.

## Usage

``` r
plot_tvcqgcomp_drift(drift_obj = NULL, observed_data = NULL, natural = NULL,
  run_obj = NULL, config = NULL, variables = NULL,
  include = c("both", "categorical", "continuous"),
  continuous_stat = c("mean", "sd", "quantile"), quantile = 0.5,
  categorical_stat = c("mean_score", "level_proportion"),
  categorical_level = NULL, ncol = NULL, observed_col = "black",
  natural_col = "grey70", observed_pch = 16, natural_lwd = 2,
  xlab = NULL, ylab = NULL, main = NULL, legend = TRUE,
  legend_pos = "bottom", ...)
```

## Arguments

- drift_obj:

  Optional object returned by
  [`diagnose_tvcqgcomp_drift()`](https://JuwelRana19.github.io/tvcQGComp/reference/diagnose_tvcqgcomp_drift.md).

- observed_data:

  Observed longitudinal data. Used when `drift_obj` is not supplied.

- natural:

  Natural-course simulated data. Used when `drift_obj` is not supplied.

- run_obj:

  Optional object returned by
  [`tvcQGComp_survival()`](https://JuwelRana19.github.io/tvcQGComp/reference/run_tvcqgcomp.md).
  Used when `drift_obj` is not supplied.

- config:

  Optional configuration object. Used when `drift_obj` is not supplied.

- variables:

  Optional character vector of variables to plot when `drift_obj` is not
  supplied.

- include:

  Whether to include `"categorical"`, `"continuous"`, or `"both"`
  variable types.

- continuous_stat:

  Continuous summary to plot: `"mean"`, `"sd"`, or `"quantile"`.

- quantile:

  Quantile to plot when `continuous_stat = "quantile"`.

- categorical_stat:

  Categorical summary to plot: `"mean_score"` or `"level_proportion"`.

- categorical_level:

  Optional category level to plot when
  `categorical_stat = "level_proportion"`. If omitted, the most common
  observed level is used for each variable.

- ncol:

  Number of columns in the panel layout. Defaults to a square-ish
  layout.

- observed_col:

  Color for observed points.

- natural_col:

  Color for natural-course predicted lines.

- observed_pch:

  Point symbol for observed values.

- natural_lwd:

  Line width for natural-course predictions.

- xlab:

  Optional x-axis label for each panel.

- ylab:

  Optional y-axis label for each panel.

- main:

  Optional overall figure title.

- legend:

  Logical; whether to add a legend.

- legend_pos:

  Legend position passed to
  [`legend`](https://rdrr.io/r/graphics/legend.html).

- ...:

  Additional arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Returns the drift object invisibly after drawing the plots.
