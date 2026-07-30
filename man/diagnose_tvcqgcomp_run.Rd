\name{diagnose_tvcqgcomp_run}
\alias{diagnose_tvcqgcomp_run}
\title{Diagnose a tvcQGcomp Run}
\description{
Run a compact set of structural and scientific checks on a \code{tvcQGcomp} result object.
The helper checks Monte Carlo panel balance, duplicate subject-time rows, finite hazards,
cumulative risk bounds, monotone cumulative risk, constant exposures within each intervention scenario,
and the final-time risk ordering across scenarios.
}
\usage{
diagnose_tvcqgcomp_run(run_obj, config = NULL, tolerance = 1e-12)
}
\arguments{
  \item{run_obj}{Object returned by \code{tvcQGComp_survival()}.}
  \item{config}{Optional configuration object. If omitted, \code{run_obj$config} is used.}
  \item{tolerance}{Numerical tolerance for monotonicity and range checks.}
}
\value{
A list with summary flags and detailed tables for Monte Carlo checks, natural course checks,
intervention checks, exposure-constancy checks, and final-time risk summaries.
}
