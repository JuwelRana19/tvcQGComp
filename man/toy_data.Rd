\name{toy_data}
\alias{toy_data}
\docType{data}
\title{Toy Longitudinal Survival Data}
\description{
A toy person-period dataset for demonstrating time-varying quantile
g-computation with a five-component exposure mixture.
}
\usage{data(toy_data)}
\format{
A \code{data.frame} with 5,584 rows and 28 variables:
\describe{
  \item{\code{UniqID}}{Simulated participant identifier.}
  \item{\code{BC}, \code{NIT}, \code{SO4}, \code{NH4}, \code{OM}}{
  Continuous time-varying mixture components in arbitrary simulated
  concentration units.}
  \item{\code{Time}, \code{TimeInn}, \code{TimeOut}}{
  Follow-up interval variables covering times 0 through 12.}
  \item{\code{status}}{Interval event indicator.}
  \item{\code{sex}, \code{visible_minority}, \code{IndigenousIdentity},
  \code{marsth}, \code{immigration_status}, \code{landed_migration},
  \code{Occupation}, \code{employment}, \code{Education}}{
  Simulated demographic and socioeconomic covariates.}
  \item{\code{age}, \code{income_inadequacy}}{
  Covariates that vary over follow-up in this toy dataset.}
  \item{\code{CanadianRegion}, \code{Urban_form}, \code{CSize},
  \code{dependency}, \code{deprivation}, \code{ethnicconcentration},
  \code{instability}}{Time-varying categorical covariates.}
}
}
\details{
The dataset contains 500 simulated individuals followed for 1 to 13 years.
Exactly 100 participants die during follow-up, corresponding to 20 percent
cumulative mortality. It is provided solely for examples and testing and
contains no real participant records. The exposure variables remain
continuous in the public dataset and are converted to quartiles internally
by \code{tvcQGComp_survival()}. Precomputed lag columns are intentionally
omitted; set \code{auto_history = TRUE} in the configuration to generate
required lag histories. An RDS copy is installed at
\code{system.file("extdata", "toy_data_500.RDS", package = "tvcQGComp")}.
}
\source{
Simulated for the tvcQGComp methods project.
}
\keyword{datasets}
