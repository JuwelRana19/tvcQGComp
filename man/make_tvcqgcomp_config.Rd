\name{make_tvcqgcomp_config}
\alias{make_tvcqgcomp_config}
\alias{generate_intervention_scenarios}
\alias{validate_tvcqgcomp_data}
\alias{make_q_scenarios}
\title{Create a tvcQGcomp Configuration and Quantized Scenarios}
\description{
Create the configuration object used by \code{tvcQGcomp}, generate simple static Q1 to QK intervention scenarios, and validate observed data against
the package requirements.
}
\usage{
make_tvcqgcomp_config(
  id,
  time_in = NULL,
  time_name = NULL,
  time_out,
  outcome,
  outcome_type = "survival",
  exposures,
  exposure_lags = NULL,
  exposure_formulas = NULL,
  exposure_types = NULL,
  time_varying_covariates,
  tvc_lags = NULL,
  covtypes = NULL,
  passive_tvc_vars = NULL,
  time_fixed_covariates,
  factor_vars = NULL,
  outcome_formula,
  tvc_formulas,
  q = NULL,
  q_levels = NULL,
  natural_course = "calibration",
  auto_history = TRUE,
  baselags = FALSE,
  meta_target = "HR",
  categorical_sim = "stochastic",
  exposome = NULL,
  exposure_scale = NULL,
  bounded_supports = NULL
)

make_q_scenarios(config, levels = NULL, prefix = "Q")

generate_intervention_scenarios(config, levels = NULL, prefix = "Q")

validate_tvcqgcomp_data(data, config)
}
\arguments{
  \item{id}{Subject identifier column name.}
  \item{time_in}{Start-of-interval time column name.}
  \item{time_name}{Alternative name for \code{time_in}.}
  \item{time_out}{End-of-interval time column name.}
  \item{outcome}{Outcome column name.}
  \item{outcome_type}{Outcome type. The current public release supports \code{"survival"} for pooled person-time survival analyses.}
  \item{exposures}{Exposure column names.}
  \item{exposure_lags}{Lagged exposure column names. Defaults to \code{paste0(exposures, "_lag1")}.}
  \item{exposure_formulas}{Optional formulas for modeled natural-course exposures.}
  \item{exposure_types}{Named vector of exposure model types. Supported values include \code{"categorical"}, \code{"multinomial"}, \code{"binary"}, \code{"normal"}, \code{"bounded_normal"}, \code{"bounded_normal_snap"}, and ordinal-logit aliases \code{"ordinal"}, \code{"ordinal_logit"}, \code{"ordinal_polr"}, \code{"ordered_logit"}, \code{"polr"}.}
  \item{time_varying_covariates}{Time-varying confounder column names.}
  \item{tvc_lags}{Lagged time-varying confounder column names.}
  \item{covtypes}{Named vector of time-varying confounder model types. Supported values include \code{"categorical"}, \code{"multinomial"}, \code{"binary"}, \code{"normal"}, \code{"bounded_normal"}, \code{"bounded_normal_snap"}, and ordinal-logit aliases \code{"ordinal"}, \code{"ordinal_logit"}, \code{"ordinal_polr"}, \code{"ordered_logit"}, \code{"polr"}.}
  \item{passive_tvc_vars}{Optional subset of \code{time_varying_covariates} to carry forward from the sampled observed id-time path without fitting a stochastic TVC model. Passive TVCs remain available to outcome, exposure, and other confounder formulas at each time point.}
  \item{time_fixed_covariates}{Baseline covariates that stay fixed within subject.}
  \item{factor_vars}{Variables that should be treated as factors during simulation and prediction.}
  \item{outcome_formula}{Outcome model formula.}
  \item{tvc_formulas}{Named list of time-varying confounder model formulas.}
  \item{q}{Primary quantization level argument. If \code{NULL}, the function falls back to \code{q_levels} (or 4).}
  \item{q_levels}{Number of quantile levels used to define the intervention grid.}
  \item{natural_course}{One or more natural-course modes. The preferred names are \code{"none"}, \code{"calibration"}, \code{"observed_exposome"}, and \code{"modeled_exposome"}. Backward-compatible aliases \code{"observed"} and \code{"modeled"} are also accepted and are normalized to \code{"calibration"} and \code{"modeled_exposome"}, respectively. \code{"calibration"} predicts outcomes on the observed data path without forward simulation; \code{"observed_exposome"} forward-simulates TVCs while keeping sampled observed exposure histories fixed; and \code{"modeled_exposome"} forward-simulates both TVCs and exposures.}
  \item{auto_history}{Logical; whether to generate required lag histories from their base variables during data preparation. When \code{TRUE}, lag columns referenced by the configuration or model formulas do not need to be present in the input data.}
  \item{baselags}{Logical scalar following the \code{gfoRmula} convention for early lag histories when pre-baseline rows are unavailable. If \code{FALSE}, early \code{*_lag2}, \code{*_lag3}, etc. values are set to 0 for numeric variables or the reference level for factors. If \code{TRUE}, early lag values are set to the subject's baseline value.}
  \item{meta_target}{Second-stage summary target. Use \code{"HR"} for the pooled hazard-based person-time multiplicative meta model, \code{"HD"} for the pooled hazard-based additive meta model, \code{"RR"} for the end-of-follow-up cumulative-risk ratio meta model, or \code{"RD"} for the end-of-follow-up cumulative-risk difference meta model. Multiple targets can be requested, for example \code{c("HR", "HD")} or \code{c("RR", "RD")}; the legacy shortcut \code{"both"} expands to \code{c("HR", "RR")}. Historical aliases \code{"hazard"} and \code{"cumrisk_final"} are also accepted.}
  \item{categorical_sim}{Categorical simulation rule for time-varying covariates and modeled exposures. \code{"stochastic"} uses multinomial probability draws; \code{"class"} uses deterministic class prediction.}
  \item{exposome}{Exposure representation. Use \code{"quantized"} for the current public release.}
  \item{exposure_scale}{Backward-compatible alias for \code{exposome}.}
  \item{bounded_supports}{Optional named list specifying legal support grids for bounded-normal snap/ordinal variables. Each entry can be a numeric vector (for example \code{1:5}), \code{list(values = ...)}, \code{list(min=, max=, by=)}, or \code{"observed"} to use observed unique values.}
  \item{data}{Observed data to validate.}
  \item{config}{A \code{tvcQGcomp} configuration object.}
  \item{levels}{Intervention levels to assign to all exposures. Defaults to \code{0:(q_levels - 1)}.}
  \item{prefix}{Scenario name prefix.}
}
\value{
\code{make_tvcqgcomp_config()} returns a named list used by the package.
\code{make_q_scenarios()} and \code{generate_intervention_scenarios()} return a named list of static intervention vectors.
\code{validate_tvcqgcomp_data()} returns \code{TRUE} invisibly and throws an error if the input data do not satisfy the required structure.
}

