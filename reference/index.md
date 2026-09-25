# Package index

## Testing for mediation

The power-enhanced global test and its outcome-specific workers.

- [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
  : Power-enhanced test for high-dimensional mediation
- [`pe_mediate()`](https://yelleknek.github.io/POEM/reference/pe_mediate.md)
  : Power-enhanced mediation test from a data frame (formula interface)
- [`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md)
  : Power-enhanced mediation test for a continuous outcome
- [`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md)
  : Power-enhanced mediation test for a binary outcome
- [`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md)
  : Power-enhanced mediation test for a count outcome

## Identifying the active mediators

- [`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md)
  : Per-mediator detail behind a power-enhanced mediation test
- [`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md)
  : Active mediators identified under each multiplicity method

## Planning and Monte Carlo studies

Size, power, identification, and sample size planning by simulation.

- [`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
  : Monte Carlo size and power curve for the PE mediation tests
- [`pe_simulation_study()`](https://yelleknek.github.io/POEM/reference/pe_simulation_study.md)
  : Run the power study under several mediation patterns at once
- [`pe_identification_study()`](https://yelleknek.github.io/POEM/reference/pe_identification_study.md)
  : Monte Carlo study of individual-mediator identification (FWER / FDR)
- [`ss_power_pe_mediation()`](https://yelleknek.github.io/POEM/reference/ss_power_pe_mediation.md)
  : Simulation-based sample size planning for the power-enhanced
  mediation test

## Simulating mediation data

- [`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md)
  : Simulate high-dimensional mediation data
- [`simulate_guo_mediation()`](https://yelleknek.github.io/POEM/reference/simulate_guo_mediation.md)
  : Simulate the real-data-motivated heterogeneous mediation setting
- [`mediation_ar1_cov()`](https://yelleknek.github.io/POEM/reference/mediation_ar1_cov.md)
  : Autoregressive (AR(1)) covariance matrix for mediator simulation

## The WHO health-expenditure analysis

- [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md)
  : Reproduce the manuscript's empirical mediation analysis
- [`WHO_mediation_design()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_design.md)
  : Build a mediation design from the WHO health-expenditure data

## Results, printing, and tidiers

- [`format(`*`<poem_tbl>`*`)`](https://yelleknek.github.io/POEM/reference/poem_tbl.md)
  [`print(`*`<poem_tbl>`*`)`](https://yelleknek.github.io/POEM/reference/poem_tbl.md)
  : Tidy printing for POEM result tables
- [`summary(`*`<poem_tbl>`*`)`](https://yelleknek.github.io/POEM/reference/summary.poem_tbl.md)
  : Summarize a POEM table
- [`plot(`*`<poem_tbl>`*`)`](https://yelleknek.github.io/POEM/reference/plot.poem_tbl.md)
  : Plot a POEM result
- [`tidy(`*`<poem_tbl>`*`)`](https://yelleknek.github.io/POEM/reference/poem_broom.md)
  [`glance(`*`<poem_tbl>`*`)`](https://yelleknek.github.io/POEM/reference/poem_broom.md)
  : Broom verbs for POEM results

## Data sets

- [`WHO_health_mediation`](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md)
  : WHO / World Bank health-expenditure mediation data (benchmark data
  set)
- [`WHO_indicator_codebook`](https://yelleknek.github.io/POEM/reference/WHO_indicator_codebook.md)
  : Codebook for the WHO health-expenditure indicators
- [`guo_calibration`](https://yelleknek.github.io/POEM/reference/guo_calibration.md)
  : Calibration constants for the real-data-motivated heterogeneous
  setting
- [`example_data`](https://yelleknek.github.io/POEM/reference/example_data.md)
  [`example_continuous`](https://yelleknek.github.io/POEM/reference/example_data.md)
  [`example_binary`](https://yelleknek.github.io/POEM/reference/example_data.md)
  [`example_count`](https://yelleknek.github.io/POEM/reference/example_data.md)
  : Example high-dimensional mediation data sets

## Package overview

- [`POEM`](https://yelleknek.github.io/POEM/reference/POEM-package.md)
  [`POEM-package`](https://yelleknek.github.io/POEM/reference/POEM-package.md)
  : POEM: Power Enhancement for Mediation in R
