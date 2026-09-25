# Plot a POEM result

A base-graphics plot for the two kinds of table POEM returns. For a
power curve from
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
it draws the empirical rejection rate of the benchmark and
power-enhanced tests against the signal-strength grid, with the nominal
level marked. For a single
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
fit it draws the per-mediator screening evidence (\\-\log\_{10}\\ of the
screening p-value) for each selected candidate mediator, with the active
ones highlighted, so it is clear which mediators drove the result.

## Usage

``` r
# S3 method for class 'poem_tbl'
plot(x, ...)
```

## Arguments

- x:

  A `poem_tbl` from
  [`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
  or
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md).

- ...:

  Further graphical parameters passed to the underlying plot.

## Value

`x`, invisibly. Called for the plot it draws.

## See also

[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md),
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md).

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEM/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEM/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEM/reference/pe_identification_study.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md),
[`pe_simulation_study()`](https://yelleknek.github.io/POEM/reference/pe_simulation_study.md),
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEM/reference/simulate_guo_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md),
[`ss_power_pe_mediation()`](https://yelleknek.github.io/POEM/reference/ss_power_pe_mediation.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                             pattern = "contrasting", c1 = 1)
fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
plot(fit)

```
