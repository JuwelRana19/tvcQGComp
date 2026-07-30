\name{plot_cumulative_risk_trajectory}
\alias{plot_cumulative_risk_trajectory}
\title{Plot Cumulative Risk Trajectories Over Time}
\description{
Draw cumulative-risk trajectories over time for natural-course and intervention scenarios from a
\code{tvcQGcomp} run object or a pre-summarized trajectory table.
This is the primary plotting function for cumulative incidence.
}
\usage{
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
}
\arguments{
  \item{x}{A \code{tvcQGcomp} run object or a data.frame/data.table containing summarized risk trajectories.}
  \item{include_natural}{Logical; when \code{x} is a run object, include the natural-course curve if available.}
  \item{natural_label}{Label to use for the natural-course curve in the legend.}
  \item{include_average}{Logical; add a line for the average of \code{average_scenarios} at each follow-up time.}
  \item{average_scenarios}{Character vector of scenario names to average. Defaults to \code{Q0} through \code{Q3}.}
  \item{average_label}{Legend label for the averaged trajectory.}
  \item{average_col}{Color for the averaged trajectory when \code{cols} is not supplied.}
  \item{scenario_col}{Name of the scenario column.}
  \item{time_col}{Name of the follow-up time column.}
  \item{risk_col}{Name of the cumulative risk column to plot.}
  \item{percent}{Logical; if \code{TRUE}, plot risk on the percent scale.}
  \item{scenario_order}{Optional character vector giving the desired plotting order of scenarios.}
  \item{scenario_labels}{Optional replacement labels for scenarios. May be a character vector named by scenario or an unnamed vector with one label per scenario.}
  \item{cols}{Optional vector of plotting colors. May be named by scenario.}
  \item{include_points}{Logical; whether to add points at each time value.}
  \item{pch}{Point character used when \code{include_points = TRUE}.}
  \item{lwd}{Line width.}
  \item{xlab}{X-axis label.}
  \item{ylab}{Y-axis label.}
  \item{main}{Plot title.}
  \item{legend_pos}{Legend position passed to \code{\link[graphics]{legend}}.}
  \item{...}{Additional arguments passed to \code{\link[graphics]{plot}}.}
}
\value{
The plotting helper returns the plotted data invisibly.
}
\examples{
\dontrun{
res <- tvcQGComp_survival(NO2_dt, cfg, mc_size = 100, seed = 1234)
plot_cumulative_risk_trajectory(res)
plot_cumulative_risk_trajectory(res, include_natural = FALSE, include_average = TRUE)
}
}
