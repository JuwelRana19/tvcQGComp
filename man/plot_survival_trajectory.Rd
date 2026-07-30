\name{plot_survival_trajectory}
\alias{plot_survival_trajectory}
\title{Plot Survival Trajectories Over Time}
\description{
Draw survival trajectories over time for natural-course and intervention
scenarios from a \code{tvcQGcomp} run object or a pre-summarized trajectory
table.
}
\usage{
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
  ylab = if (percent) "Survival (\%)" else "Survival probability",
  main = "Survival over time",
  legend_pos = "bottomleft",
  ...
)
}
\arguments{
  \item{x}{A \code{tvcQGcomp} run object or a data frame containing summarized trajectories.}
  \item{include_natural}{Logical; include the natural-course curve when available.}
  \item{natural_label}{Legend label for the natural-course curve.}
  \item{include_average}{Logical; add a line for the average of \code{average_scenarios} at each follow-up time.}
  \item{average_scenarios}{Character vector of scenario names to average. Defaults to \code{Q0} through \code{Q3}.}
  \item{average_label}{Legend label for the averaged trajectory.}
  \item{average_col}{Color for the averaged trajectory when \code{cols} is not supplied.}
  \item{scenario_col}{Name of the scenario column.}
  \item{time_col}{Name of the follow-up time column.}
  \item{survival_col}{Name of the survival column to plot.}
  \item{percent}{Logical; plot survival on the percentage scale.}
  \item{scenario_order}{Optional character vector giving scenario order.}
  \item{scenario_labels}{Optional replacement labels for scenarios.}
  \item{cols}{Optional plotting colors, optionally named by scenario.}
  \item{include_points}{Logical; add points at each follow-up time.}
  \item{pch}{Point character.}
  \item{lwd}{Line width.}
  \item{xlab}{Horizontal-axis label.}
  \item{ylab}{Vertical-axis label.}
  \item{main}{Plot title.}
  \item{legend_pos}{Legend position passed to \code{\link[graphics]{legend}}.}
  \item{...}{Additional arguments passed to \code{\link[graphics]{plot}}.}
}
\value{
Returns the plotted trajectory data invisibly.
}
\examples{
\dontrun{
fit <- tvcQGComp_survival(dat, cfg, intervention_scenarios = scenarios)
plot_survival_trajectory(fit, include_natural = FALSE)
plot_survival_trajectory(fit, include_natural = FALSE, include_average = TRUE)
}
}
