# Autoregressive (AR(1)) covariance matrix for mediator simulation

Builds the \\p \times p\\ first-order autoregressive covariance matrix
\\\Sigma = (\rho^{\|i-j\|})\_{i,j}\\, the intercorrelation structure
imposed on the mediator noise throughout the manuscript's Monte Carlo
studies. It is exported as a small reusable utility so a user can
inspect or reuse the exact covariance the simulations use.

## Usage

``` r
mediation_ar1_cov(p, rho = 0.5)
```

## Arguments

- p:

  Number of mediators (matrix dimension); a positive integer.

- rho:

  Autocorrelation parameter in \\(-1, 1)\\. The covariance between
  mediators \\i\\ and \\j\\ is \\\rho^{\|i-j\|}\\, so nearby mediators
  are more strongly correlated and the correlation decays geometrically
  with separation. Default 0.5, the value used in the manuscript.

## Value

A \\p \times p\\ numeric matrix with ones on the diagonal and
\\\rho^{\|i-j\|}\\ off the diagonal.

## See also

[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md),
which uses this covariance to draw correlated mediator noise.

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEM/reference/guo_calibration.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEM/reference/pe_identification_study.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md),
[`pe_simulation_study()`](https://yelleknek.github.io/POEM/reference/pe_simulation_study.md),
[`plot.poem_tbl()`](https://yelleknek.github.io/POEM/reference/plot.poem_tbl.md),
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEM/reference/simulate_guo_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md),
[`ss_power_pe_mediation()`](https://yelleknek.github.io/POEM/reference/ss_power_pe_mediation.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
mediation_ar1_cov(4, rho = 0.5)
#>       [,1] [,2] [,3]  [,4]
#> [1,] 1.000 0.50 0.25 0.125
#> [2,] 0.500 1.00 0.50 0.250
#> [3,] 0.250 0.50 1.00 0.500
#> [4,] 0.125 0.25 0.50 1.000
```
