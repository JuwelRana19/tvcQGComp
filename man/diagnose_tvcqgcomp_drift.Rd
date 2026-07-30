\name{diagnose_tvcqgcomp_drift}
\alias{diagnose_tvcqgcomp_drift}
\title{Diagnose Drift Between Observed and Natural-Course TVCs}
\description{
Compare observed and simulated natural-course trajectories over time for time-varying covariates.
Categorical variables are summarized through time-specific level proportions, while continuous variables
are summarized through time-specific means, standard deviations, and selected quantiles.
}
\usage{
diagnose_tvcqgcomp_drift(observed_data = NULL, natural = NULL, run_obj = NULL,
  config = NULL, variables = NULL, include = c("both", "categorical",
  "continuous"), continuous_probs = c(0.1, 0.5, 0.9))
}
\arguments{
  \item{observed_data}{Observed longitudinal data used as the calibration reference.}
  \item{natural}{Simulated natural-course data to compare against the observed data.}
  \item{run_obj}{Optional object returned by \code{tvcQGComp_survival()}. When supplied, \code{config},
  \code{observed_data}, and \code{natural} are filled from the run object when available.}
  \item{config}{Optional configuration object. If omitted, \code{run_obj$config} is used.}
  \item{variables}{Optional character vector of variables to check. Defaults to the configured time-varying covariates.}
  \item{include}{Whether to return checks for \code{"categorical"}, \code{"continuous"}, or \code{"both"} types.}
  \item{continuous_probs}{Quantiles to compare over time for continuous variables.}
}
\value{
A list with the inferred variable types and separate detailed and summary tables for categorical and
continuous drift over time.
}
