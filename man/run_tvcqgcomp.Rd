\name{run_tvcqgcomp}
\alias{tvcQGComp_survival}
\title{Run the tvcQGcomp Pipeline}
\description{
Fit the nuisance models, create the Monte Carlo skeleton, simulate the natural course and static quantized interventions,
and fit the pooled person-time meta model. The canonical user-facing survival entry point is
\code{tvcQGComp_survival()}.
}
\usage{
tvcQGComp_survival(
  data,
  config,
  mc_size = 10000L,
  intervention_scenarios = NULL,
  natural_course = config$natural_course,
  seed = 1234L,
  replace_mc = NULL,
  meta_formula = NULL,
  meta_family = NULL,
  meta_target = config$meta_target,
  fit_meta = NULL,
  fit_exposure_models = NULL,
  exposure_scale = config$exposure_scale,
  quantization_breaks = NULL,
  joint_effect_fn = NULL,
  verbose = FALSE
)
}
\arguments{
  \item{data}{Observed input data.}
  \item{config}{A \code{tvcQGcomp} configuration object.}
  \item{mc_size}{Monte Carlo sample size.}
  \item{intervention_scenarios}{Named list of static quantized intervention values.}
  \item{natural_course}{Any subset of the preferred labels \code{"calibration"}, \code{"observed_exposome"}, and \code{"modeled_exposome"}, or \code{"none"}. Backward-compatible aliases \code{"observed"} and \code{"modeled"} are also accepted. \code{"calibration"} predicts outcomes on the observed data path, \code{"observed_exposome"} forward-simulates TVCs while keeping sampled observed exposure histories, and \code{"modeled_exposome"} forward-simulates both TVCs and exposures.}
  \item{seed}{Random seed for Monte Carlo resampling and simulation.}
  \item{replace_mc}{Whether to sample subjects with replacement when creating the Monte Carlo skeleton.}
  \item{meta_formula}{Optional pooled meta-model formula.}
  \item{meta_family}{Optional model family for the pooled meta model. If a single target is requested and \code{meta_family} is \code{NULL}, defaults are \code{quasibinomial("logit")} for \code{"HR"}, \code{gaussian("identity")} for \code{"HD"}, \code{quasibinomial("identity")} for \code{"RD"}, and \code{quasibinomial("log")} for \code{"RR"}. When multiple meta targets are requested, pass either \code{NULL} to use defaults for each target or a named list such as \code{list(HR = quasibinomial("logit"), HD = gaussian("identity"))}.}
  \item{meta_target}{Second-stage summary target. Use \code{"HR"}, \code{"HD"}, \code{"RR"}, \code{"RD"}, or a vector such as \code{c("HR", "HD")}. Historical aliases \code{"hazard"} and \code{"cumrisk_final"} are also accepted; the legacy shortcut \code{"both"} expands to \code{c("HR", "RR")}. When multiple targets are requested, the first-stage simulation is run once and all requested second-stage meta-models are fit from the same simulated intervention data.}
  \item{fit_meta}{Optional logical controlling whether the second-stage meta-model is fit. Defaults to \code{TRUE} for quantized exposures. When \code{FALSE}, \code{risk_trajectory} still reports time-specific scenario risks and survival, but no HR, HD, RR, or RD effect estimate is produced.}
  \item{fit_exposure_models}{Optional logical controlling whether exposure evolution models are fit. Defaults to \code{FALSE} when \code{natural_course = "none"} and every intervention scenario deterministically assigns every exposure; otherwise defaults to \code{TRUE}.}
  \item{exposure_scale}{Exposure representation used during first-stage fitting and simulation. For this minimal release, use \code{"quantized"}.}
  \item{quantization_breaks}{Optional fixed cutpoints used when \code{exposure_scale = "quantized"}.}
  \item{joint_effect_fn}{Optional function mapping a scenario vector to a scalar joint effect.}
  \item{verbose}{Logical; print progress and status messages during fitting and simulation.}
}
\value{
A list containing \code{mc_data}, \code{natural}, \code{interventions}, \code{intervention_data}, \code{risk_trajectory}, \code{meta_model}, \code{meta_models}, \code{meta_effect_summary}, calibration outputs, and the fitted nuisance models. \code{risk_trajectory} is the sole public scenario-prediction table and reports time-specific cumulative risk and survival by intervention scenario. HR, HD, RR, and RD estimates are derived exclusively from the requested second-stage meta-models and reported in \code{meta_effect_summary}, which contains only \code{meta_target}, \code{estimate}, and \code{increment}. For HR and RR, \code{estimate} is the exponentiated \code{joint_effect} coefficient; for HD and RD, it is the untransformed coefficient. Each estimate represents a one-unit increase in the joint quantile intervention. When \code{fit_meta = FALSE}, the meta-model elements and \code{meta_effect_summary} are \code{NULL}. When \code{natural_course = "none"}, \code{natural} and the natural-course calibration outputs are \code{NULL}. The Monte Carlo skeleton inherits baseline and lagged predictors from the observed data, so simulated baseline/history support must be complete for the predictors required by the nuisance models. The \code{calibration} element includes \code{calibration}/\code{calibration_curve}, \code{observed_exposome}/\code{observed_exposome_curve}, and \code{modeled_exposome}/\code{modeled_exposome_curve}; backward-compatible aliases \code{observed} and \code{modeled} are retained.
}

