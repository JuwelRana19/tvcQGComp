\name{as_risk_curve_df}
\alias{as_risk_curve_df}
\title{Extract Tidy Risk-Curve Summaries}
\description{
Combine natural-course and intervention risk-trajectory summaries from a \code{tvcQGcomp}
run object into one tidy table for downstream plotting or tabulation.
}
\usage{
as_risk_curve_df(
  run_obj,
  include_natural = TRUE,
  natural_label = "Natural Course",
  time_col = NULL,
  scenario_col = "scenario"
)
}
\arguments{
  \item{run_obj}{A \code{tvcQGcomp} run object containing \code{natural_curve}
  and/or \code{risk_trajectory}. Older objects using \code{intervention_curve}
  remain supported.}
  \item{include_natural}{Logical; include the natural-course summary if available.}
  \item{natural_label}{Scenario label to assign to the natural-course rows.}
  \item{time_col}{Optional time column name. Defaults to \code{run_obj$config$time_in}.}
  \item{scenario_col}{Name of the scenario column in the returned table.}
}
\value{
A \code{data.table} containing one row per scenario and time point, with columns such as
\code{mean_risk}, \code{mean_risk_per_1000}, \code{mean_survival}, and
\code{mean_survival_per_1000}.
}
\details{
This helper is useful when you want to reproduce a \code{gfoRmula}-style risk-curve plot
using \code{ggplot2} or custom tabulation code.
}
\examples{
\dontrun{
curve_dt <- as_risk_curve_df(res)

ggplot(curve_dt, aes(x = TimeInn, y = mean_risk, color = scenario)) +
  geom_line()
}
}
