\name{bootstrap_tvcqgcomp}
\alias{tvcQGComp_survival_boot}
\alias{summarize_tvcqgcomp_bootstrap}
\title{Bootstrap tvcQGcomp Fits}
\description{
Run a subject-level cluster bootstrap for \code{tvcQGcomp}. In regular mode the function returns
replicate-level meta-model coefficients. In audit mode it also stores the per-replicate bootstrap data,
Monte Carlo data, natural course output, intervention outputs, and risk curves. The canonical user-facing
survival bootstrap entry point is \code{tvcQGComp_survival_boot()}.
}
\usage{
tvcQGComp_survival_boot(
  data,
  config,
  n_boot = 200L,
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
  parallel = FALSE,
  n_workers = 2L,
  batch_size = NULL,
  checkpoint_file = NULL,
  resume_from_checkpoint = TRUE,
  verbose = TRUE,
  stop_on_error = TRUE,
  audit_mode = FALSE,
  keep_boot_objects = FALSE,
  keep_boot_curves = FALSE
)

summarize_tvcqgcomp_bootstrap(boot_res)
}
\arguments{
  \item{data, config, mc_size, intervention_scenarios, natural_course, seed, replace_mc, meta_formula, meta_family, meta_target, fit_meta, fit_exposure_models, exposure_scale, quantization_breaks, joint_effect_fn}{See \code{tvcQGComp_survival()}. When multiple meta targets are requested, bootstrap result columns are prefixed by target, for example \code{coef_HR_joint_effect}, \code{coef_HD_joint_effect}, and \code{coef_RR_joint_effect}. When \code{fit_meta = FALSE}, replicate outputs contain scenario-specific final risks instead of meta-model coefficients.}
  \item{n_boot}{Number of bootstrap replicates.}
  \item{parallel}{Logical; whether to evaluate bootstrap replicates in parallel.}
  \item{n_workers}{Number of workers when \code{parallel = TRUE}.}
  \item{batch_size}{Number of replicates per parallel batch.}
  \item{checkpoint_file}{Optional RDS file to update during the bootstrap run.}
  \item{resume_from_checkpoint}{Logical; if \code{TRUE} and \code{checkpoint_file} exists, resume from the last completed bootstrap replicate stored in the checkpoint.}
  \item{verbose}{Logical; print bootstrap progress.}
  \item{stop_on_error}{Logical; if \code{TRUE}, stop the bootstrap at the first failed replicate after saving a checkpoint. If \code{FALSE}, failed replicates are retained with an \code{error: ...} status and the bootstrap continues.}
  \item{audit_mode}{Logical; if \code{TRUE}, keep replicate objects for auditing.}
  \item{keep_boot_objects}{Logical; explicit switch to retain audit objects.}
  \item{keep_boot_curves}{Logical; retain a lightweight per-replicate curve table containing time-specific risk and survival without storing the full bootstrap audit objects.}
  \item{boot_res}{Result from \code{tvcQGComp_survival_boot()}.}
}
\value{
\code{tvcQGComp_survival_boot()} returns either a data table of replicate meta-model coefficients or scenario-specific final risks or,
when curve summaries and/or audit objects are retained, a list with \code{results} and optional \code{curve_results},
\code{audit_objects}, and \code{audit_mode}.
\code{summarize_tvcqgcomp_bootstrap()} returns a summary table with means,
standard errors, and percentile intervals.
}

