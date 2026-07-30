# Plot Cumulative-Risk Trajectories

Plot cumulative incidence from the time-specific intervention
predictions stored in a fitted `tvcQGComp` run object.

## Usage

``` r
# S3 method for class 'tvcqgcomp_run'
plot(x, ...)
```

## Arguments

- x:

  A fitted `tvcqgcomp_run` object.

- ...:

  Arguments passed to
  [`plot_cumulative_risk_trajectory()`](https://JuwelRana19.github.io/tvcQGComp/reference/plot_cumulative_risk_trajectory.md).

## Details

`plot(x)` is a convenience alias for
`plot_cumulative_risk_trajectory(x)`. Use `plot_survival_trajectory(x)`
when survival curves are required.

## Value

Returns the plotted trajectory data invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
plot(fit, include_natural = FALSE)
plot_survival_trajectory(fit, include_natural = FALSE)
} # }
```
