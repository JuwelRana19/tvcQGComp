# Package index

## Configure and validate

Define the longitudinal data structure and quantized interventions.

- [`make_tvcqgcomp_config()`](https://juwelrana.com/tvcQGComp/reference/make_tvcqgcomp_config.md)
  [`make_q_scenarios()`](https://juwelrana.com/tvcQGComp/reference/make_tvcqgcomp_config.md)
  [`generate_intervention_scenarios()`](https://juwelrana.com/tvcQGComp/reference/make_tvcqgcomp_config.md)
  [`validate_tvcqgcomp_data()`](https://juwelrana.com/tvcQGComp/reference/make_tvcqgcomp_config.md)
  : Create a tvcQGcomp Configuration and Quantized Scenarios
- [`toy_data`](https://juwelrana.com/tvcQGComp/reference/toy_data.md) :
  Toy Longitudinal Survival Data

## Estimate effects

Fit the main survival analysis and obtain bootstrap uncertainty.

- [`tvcQGComp_survival()`](https://juwelrana.com/tvcQGComp/reference/run_tvcqgcomp.md)
  : Run the tvcQGcomp Pipeline
- [`tvcQGComp_survival_boot()`](https://juwelrana.com/tvcQGComp/reference/bootstrap_tvcqgcomp.md)
  [`summarize_tvcqgcomp_bootstrap()`](https://juwelrana.com/tvcQGComp/reference/bootstrap_tvcqgcomp.md)
  : Bootstrap tvcQGcomp Fits

## Summarize and plot

Extract and visualize risk trajectories.

- [`as_risk_curve_df()`](https://juwelrana.com/tvcQGComp/reference/as_risk_curve_df.md)
  : Extract Tidy Risk-Curve Summaries
- [`plot(`*`<tvcqgcomp_run>`*`)`](https://juwelrana.com/tvcQGComp/reference/plot.tvcqgcomp_run.md)
  : Plot Cumulative-Risk Trajectories
- [`plot_cumulative_risk_trajectory()`](https://juwelrana.com/tvcQGComp/reference/plot_cumulative_risk_trajectory.md)
  : Plot Cumulative Risk Trajectories Over Time
- [`plot_survival_trajectory()`](https://juwelrana.com/tvcQGComp/reference/plot_survival_trajectory.md)
  : Plot Survival Trajectories Over Time

## Diagnose models

Check fitted runs and natural-course trajectory drift.

- [`diagnose_tvcqgcomp_run()`](https://juwelrana.com/tvcQGComp/reference/diagnose_tvcqgcomp_run.md)
  : Diagnose a tvcQGcomp Run
- [`diagnose_tvcqgcomp_drift()`](https://juwelrana.com/tvcQGComp/reference/diagnose_tvcqgcomp_drift.md)
  : Diagnose Drift Between Observed and Natural-Course TVCs
- [`plot_tvcqgcomp_drift()`](https://juwelrana.com/tvcQGComp/reference/plot_tvcqgcomp_drift.md)
  : Plot Observed vs Natural-Course Drift Over Time
