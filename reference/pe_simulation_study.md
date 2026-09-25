# Run the power study under several mediation patterns at once

A convenience wrapper around
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
that runs the size and power study under more than one mediation pattern
in a single call and stacks the results, so a homogeneous and a
contrasting curve (the two panels of a figure in the article come back
together. Each pattern is swept over the same signal grid.

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
  alpha_level = 0.05,
  method = c("Bonferroni", "BH", "BY"),
  error_level = 0.05,
  lambda_grid = NULL,
  lambda_grid_reduced = NULL,
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

- c1_grid, c2, n_rep, alpha_level, method, error_level, lambda_grid,
  lambda_grid_reduced, cores, seed:

  Passed to
  [`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md).

## Value

A tidy `data.frame` (class `poem_tbl`) with a leading `pattern` column
and the
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
columns (`c1`, `rejection_hdmm`, `rejection_pe`, `n_valid`, `n_empty`)
for each pattern. The settings are recorded in the attributes `outcome`,
`n`, `p`, `n_rep`, `alpha_level`, `error_level`, and `lambda_grid`.

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
# n_rep = 5 keeps this example fast; a rate from five replications has a
# Monte Carlo standard error of up to 0.22, so read the shape, not the
# numbers. A reported study uses the article's design (n = 300, p = 500)
# and 1000 replications.
set.seed(113)
pe_simulation_study(n = 120, p = 50, outcome = "continuous",
                    c1_grid = c(0, 1), n_rep = 5)
#>  pattern     c1 rejection_hdmm rejection_pe n_valid n_empty
#>  homogeneous 0  0              0            5       0      
#>  homogeneous 1  0.6            1            5       0      
#>  contrasting 0  0              0            5       0      
#>  contrasting 1  0              0.4          5       0      
#> 
#> Outcome: continuous
```
