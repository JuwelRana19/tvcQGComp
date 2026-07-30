# Simulated Longitudinal Survival Data

A simulated person-period dataset for demonstrating time-varying
quantile g-computation with a five-component exposure mixture.

## Usage

``` r
data(sim_data)
```

## Format

A `data.frame` with 5,584 rows and 28 variables:

- `UniqID`:

  Simulated participant identifier.

- `BC`, `NIT`, `SO4`, `NH4`, `OM`:

  Continuous time-varying mixture components in arbitrary simulated
  concentration units.

- `Time`, `TimeInn`, `TimeOut`:

  Follow-up interval variables covering times 0 through 12.

- `status`:

  Interval event indicator.

- `sex`, `visible_minority`, `IndigenousIdentity`, `marsth`,
  `immigration_status`, `landed_migration`, `Occupation`, `employment`,
  `Education`:

  Simulated demographic and socioeconomic covariates.

- `age`, `income_inadequacy`:

  Covariates that vary over follow-up in this simulated dataset.

- `CanadianRegion`, `Urban_form`, `CSize`, `dependency`, `deprivation`,
  `ethnicconcentration`, `instability`:

  Time-varying categorical covariates.

## Details

The dataset contains 500 simulated individuals followed for 1 to 13
years. Exactly 100 participants die during follow-up, corresponding to
20 percent cumulative mortality. It is provided solely for examples and
testing and contains no real participant records. The exposure variables
remain continuous in the public dataset and are converted to quartiles
internally by
[`tvcQGComp_survival()`](https://JuwelRana19.github.io/tvcQGComp/reference/run_tvcqgcomp.md).
Precomputed lag columns are intentionally omitted; set
`auto_history = TRUE` in the configuration to generate required lag
histories. An RDS copy is installed at
`system.file("extdata", "sim_data_500.RDS", package = "tvcQGComp")`.

## Source

Simulated for the tvcQGComp methods project.
