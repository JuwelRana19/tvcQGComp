# Diagnose a tvcQGcomp Run

Run a compact set of structural and scientific checks on a `tvcQGcomp`
result object. The helper checks Monte Carlo panel balance, duplicate
subject-time rows, finite hazards, cumulative risk bounds, monotone
cumulative risk, constant exposures within each intervention scenario,
and the final-time risk ordering across scenarios.

## Usage

``` r
diagnose_tvcqgcomp_run(run_obj, config = NULL, tolerance = 1e-12)
```

## Arguments

- run_obj:

  Object returned by
  [`tvcQGComp_survival()`](https://JuwelRana19.github.io/tvcQGComp/reference/run_tvcqgcomp.md).

- config:

  Optional configuration object. If omitted, `run_obj$config` is used.

- tolerance:

  Numerical tolerance for monotonicity and range checks.

## Value

A list with summary flags and detailed tables for Monte Carlo checks,
natural course checks, intervention checks, exposure-constancy checks,
and final-time risk summaries.
