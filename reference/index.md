# Package index

## Configure and validate

Define the longitudinal data structure and quantized interventions.

- [`make_tvcqgcomp_config()`](https://JuwelRana19.github.io/tvcQGComp/reference/make_tvcqgcomp_config.md)
  [`make_q_scenarios()`](https://JuwelRana19.github.io/tvcQGComp/reference/make_tvcqgcomp_config.md)
  [`generate_intervention_scenarios()`](https://JuwelRana19.github.io/tvcQGComp/reference/make_tvcqgcomp_config.md)
  [`validate_tvcqgcomp_data()`](https://JuwelRana19.github.io/tvcQGComp/reference/make_tvcqgcomp_config.md)
  : Create a tvcQGcomp Configuration and Quantized Scenarios
- [`toy_data`](https://JuwelRana19.github.io/tvcQGComp/reference/toy_data.md)
  : Toy Longitudinal Survival Data

## Estimate effects

Fit the main survival analysis and obtain bootstrap uncertainty.

- [`tvcQGComp_survival()`](https://JuwelRana19.github.io/tvcQGComp/reference/run_tvcqgcomp.md)
  : Run the tvcQGcomp Pipeline
- [`tvcQGComp_survival_boot()`](https://JuwelRana19.github.io/tvcQGComp/reference/bootstrap_tvcqgcomp.md)
  [`summarize_tvcqgcomp_bootstrap()`](https://JuwelRana19.github.io/tvcQGComp/reference/bootstrap_tvcqgcomp.md)
  : Bootstrap tvcQGcomp Fits

## Summarize and plot

Extract and visualize risk trajectories.

- [`as_risk_curve_df()`](https://JuwelRana19.github.io/tvcQGComp/reference/as_risk_curve_df.md)
  : Extract Tidy Risk-Curve Summaries
- [`plot(`*`<tvcqgcomp_run>`*`)`](https://JuwelRana19.github.io/tvcQGComp/reference/plot.tvcqgcomp_run.md)
  : Plot Cumulative-Risk Trajectories
- [`plot_cumulative_risk_trajectory()`](https://JuwelRana19.github.io/tvcQGComp/reference/plot_cumulative_risk_trajectory.md)
  : Plot Cumulative Risk Trajectories Over Time
- [`plot_survival_trajectory()`](https://JuwelRana19.github.io/tvcQGComp/reference/plot_survival_trajectory.md)
  : Plot Survival Trajectories Over Time

## Diagnose models

Check fitted runs and natural-course trajectory drift.

- [`diagnose_tvcqgcomp_run()`](https://JuwelRana19.github.io/tvcQGComp/reference/diagnose_tvcqgcomp_run.md)
  : Diagnose a tvcQGcomp Run
- [`diagnose_tvcqgcomp_drift()`](https://JuwelRana19.github.io/tvcQGComp/reference/diagnose_tvcqgcomp_drift.md)
  : Diagnose Drift Between Observed and Natural-Course TVCs
- [`plot_tvcqgcomp_drift()`](https://JuwelRana19.github.io/tvcQGComp/reference/plot_tvcqgcomp_drift.md)
  : Plot Observed vs Natural-Course Drift Over Time
