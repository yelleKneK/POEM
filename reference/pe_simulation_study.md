# Run the power study under several mediation patterns at once

A convenience wrapper around
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
that runs the size and power study under more than one mediation pattern
in a single call and stacks the results, so a homogeneous and a
contrasting curve (the two panels of a figure in the manuscript come
back together. Each pattern is swept over the same signal grid.

## Usage

``` r
pe_simulation_study(
  n,
  p,
  outcome = c("continuous", "binary", "count"),
  patterns = c("homogeneous", "contrasting"),
  c1_grid = seq(0, 1, by = 0.25),
  c2 = 0.5,
  n_rep = 100,
  alpha = 0.05,
  method = c("Bonferroni", "BH", "BY"),
  error_level = 0.05,
  cores = 1L,
  seed = NULL
)
```

## Arguments

- n, p:

  Observations and candidate mediators per simulated data set.

- outcome:

  Outcome type: `"continuous"`, `"binary"`, or `"count"`.

- patterns:

  Character vector of mediation patterns to run. Default
  `c("homogeneous", "contrasting")`.

- c1_grid, c2, n_rep, alpha, method, error_level, cores, seed:

  Passed to
  [`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md).

## Value

A tidy `data.frame` (class `poem_tbl`) with a leading `pattern` column
and the
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
columns (`c1`, `rejection_hdmm`, `rejection_pe`, `n_valid`) for each
pattern.

## See also

[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
for a single pattern,
[`plot.poem_tbl()`](https://yelleknek.github.io/POEM/reference/plot.poem_tbl.md)
to draw one pattern's curves.

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEM/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEM/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEM/reference/pe_identification_study.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md),
[`plot.poem_tbl()`](https://yelleknek.github.io/POEM/reference/plot.poem_tbl.md),
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEM/reference/simulate_guo_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md),
[`ss_power_pe_mediation()`](https://yelleknek.github.io/POEM/reference/ss_power_pe_mediation.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# \donttest{
set.seed(113)
pe_simulation_study(n = 120, p = 50, outcome = "continuous",
                    c1_grid = c(0, 0.5, 1), n_rep = 20)
#>  pattern     c1  rejection_hdmm rejection_pe n_valid
#>  homogeneous 0   0.05           0.05         20     
#>  homogeneous 0.5 0.35           0.75         20     
#>  homogeneous 1   0.9            1            20     
#>  contrasting 0   0.05           0.05         20     
#>  contrasting 0.5 0.1            0.45         20     
#>  contrasting 1   0.05           0.8          20     
#> 
#> Outcome model: continuous
# }
```
