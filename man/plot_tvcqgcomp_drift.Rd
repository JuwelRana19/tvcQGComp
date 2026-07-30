\name{plot_tvcqgcomp_drift}
\alias{plot_tvcqgcomp_drift}
\title{Plot Observed vs Natural-Course Drift Over Time}
\description{
Create small-multiple calibration plots that compare observed longitudinal trajectories with
predicted natural-course trajectories from \code{tvcQGcomp}. Continuous variables can be shown
 using means, standard deviations, or selected quantiles. Categorical variables can be shown using
 either time-specific mean scores or the proportion in a selected level.
}
\usage{
plot_tvcqgcomp_drift(drift_obj = NULL, observed_data = NULL, natural = NULL,
  run_obj = NULL, config = NULL, variables = NULL,
  include = c("both", "categorical", "continuous"),
  continuous_stat = c("mean", "sd", "quantile"), quantile = 0.5,
  categorical_stat = c("mean_score", "level_proportion"),
  categorical_level = NULL, ncol = NULL, observed_col = "black",
  natural_col = "grey70", observed_pch = 16, natural_lwd = 2,
  xlab = NULL, ylab = NULL, main = NULL, legend = TRUE,
  legend_pos = "bottom", ...)
}
\arguments{
  \item{drift_obj}{Optional object returned by \code{diagnose_tvcqgcomp_drift()}.}
  \item{observed_data}{Observed longitudinal data. Used when \code{drift_obj} is not supplied.}
  \item{natural}{Natural-course simulated data. Used when \code{drift_obj} is not supplied.}
  \item{run_obj}{Optional object returned by \code{tvcQGComp_survival()}. Used when \code{drift_obj} is not supplied.}
  \item{config}{Optional configuration object. Used when \code{drift_obj} is not supplied.}
  \item{variables}{Optional character vector of variables to plot when \code{drift_obj} is not supplied.}
  \item{include}{Whether to include \code{"categorical"}, \code{"continuous"}, or \code{"both"} variable types.}
  \item{continuous_stat}{Continuous summary to plot: \code{"mean"}, \code{"sd"}, or \code{"quantile"}.}
  \item{quantile}{Quantile to plot when \code{continuous_stat = "quantile"}.}
  \item{categorical_stat}{Categorical summary to plot: \code{"mean_score"} or \code{"level_proportion"}.}
  \item{categorical_level}{Optional category level to plot when \code{categorical_stat = "level_proportion"}. If omitted, the most common observed level is used for each variable.}
  \item{ncol}{Number of columns in the panel layout. Defaults to a square-ish layout.}
  \item{observed_col}{Color for observed points.}
  \item{natural_col}{Color for natural-course predicted lines.}
  \item{observed_pch}{Point symbol for observed values.}
  \item{natural_lwd}{Line width for natural-course predictions.}
  \item{xlab}{Optional x-axis label for each panel.}
  \item{ylab}{Optional y-axis label for each panel.}
  \item{main}{Optional overall figure title.}
  \item{legend}{Logical; whether to add a legend.}
  \item{legend_pos}{Legend position passed to \code{\link[graphics]{legend}}.}
  \item{...}{Additional arguments passed to \code{\link[graphics]{plot}}.}
}
\value{
Returns the drift object invisibly after drawing the plots.
}
