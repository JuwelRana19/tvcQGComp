\name{plot.tvcqgcomp_run}
\alias{plot.tvcqgcomp_run}
\title{Plot Cumulative-Risk Trajectories}
\description{
Plot cumulative incidence from the time-specific intervention predictions
stored in a fitted \code{tvcQGComp} run object.
}
\usage{
\method{plot}{tvcqgcomp_run}(x, ...)
}
\arguments{
  \item{x}{A fitted \code{tvcqgcomp_run} object.}
  \item{...}{Arguments passed to \code{plot_cumulative_risk_trajectory()}.}
}
\details{
\code{plot(x)} is a convenience alias for
\code{plot_cumulative_risk_trajectory(x)}. Use
\code{plot_survival_trajectory(x)} when survival curves are required.
}
\value{
Returns the plotted trajectory data invisibly.
}
\examples{
\dontrun{
plot(fit, include_natural = FALSE)
plot_survival_trajectory(fit, include_natural = FALSE)
}
}
