# Simulation-based sample size planning for the power-enhanced mediation test

Estimates, by Monte Carlo simulation, the sample size needed for the
power-enhanced test to reach a target power under a given mediation
pattern and signal strength. It evaluates the empirical power of the PE
test (and, for comparison, the benchmark Wald test) at each candidate
sample size on a grid, and reports the smallest grid value that reaches
the target. This is the design counterpart of the analysis function
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md):
use it to plan a study, then
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
to analyze it.

## Usage

``` r
ss_power_pe_mediation(
  outcome = c("continuous", "binary", "count"),
  pattern = c("homogeneous", "contrasting", "heterogeneous"),
  p,
  n_grid,
  c1 = 1,
  c2 = 0.5,
  target_power = 0.8,
  n_rep = 200,
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

- outcome:

  Outcome type: `"continuous"`, `"binary"`, or `"count"`.

- pattern:

  Mediation pattern: `"homogeneous"` or `"contrasting"`
  (`"heterogeneous"` is a synonym for `"contrasting"`).

- p:

  Number of candidate mediators.

- n_grid:

  Vector of candidate sample sizes to evaluate.

- c1:

  Signal strength (the exposure-on-mediator scale). Default 1.

- c2:

  Direct effect. Default 0.5.

- target_power:

  Desired power for the PE test. Default 0.8.

- n_rep:

  Monte Carlo replications per candidate `n`. Default 200.

- alpha_level:

  Significance level. Default 0.05.

- method, error_level:

  Passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md).

- lambda_grid, lambda_grid_reduced:

  Passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md);
  the default `NULL` uses the family's default grid at each candidate
  `n` (for a continuous outcome the article's grid rescaled by
  [`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md),
  so the grid follows the sample size being planned).

- cores:

  Number of CPU cores (Unix forking). A seeded parallel run is
  reproducible across runs at the same `cores` but need not match a
  serial run. Default 1.

- seed:

  Optional integer seed, set locally and restored on exit.

## Value

A tidy `data.frame` (class `poem_tbl`) with one row per candidate sample
size and columns `n`, `power_pe`, `power_hdmm`, `n_valid`, and `n_empty`
(see
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)).
The smallest `n` reaching `target_power` for the PE test is stored in
the `"recommended_n"` attribute (`NA` if no grid value reaches it), with
`"target_power"`, `"outcome"`, `"pattern"`, `"c1"`, `"c2"`, `"n_rep"`,
and `"alpha_level"`.

## Details

This is a simulation-based procedure, not a closed-form power formula:
no analytic power expression exists for the PE test, so sample size is
instead determined by direct Monte Carlo simulation. Power is estimated
by simulating `n_rep` data sets at each candidate `n` with
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md)
and recording how often the test rejects at level `alpha_level`. Because
the estimate is a Monte Carlo proportion, it carries simulation error of
roughly \\\sqrt{power(1 - power) / n\\rep}\\; raise `n_rep` for a
smoother curve and a more stable recommendation. The grid approach
(rather than a root search) is deliberate: the power curve is monotone
but noisy, so a search can stop early on a lucky draw, whereas the grid
shows the whole trajectory.

## See also

[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md),
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md).

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEM/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEM/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEM/reference/pe_identification_study.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md),
[`pe_simulation_study()`](https://yelleknek.github.io/POEM/reference/pe_simulation_study.md),
[`plot.poem_tbl()`](https://yelleknek.github.io/POEM/reference/plot.poem_tbl.md),
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEM/reference/simulate_guo_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# n_rep = 5 keeps this example fast; a power estimate from five
# replications has a Monte Carlo standard error of up to 0.22, so the
# recommended n here is a demonstration, not a plan. Planning a study
# deserves n_rep of several hundred.
set.seed(113)
plan <- ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
                              p = 50, n_grid = c(80, 120, 160), n_rep = 5)
plan
#>  n   power_pe power_hdmm n_valid n_empty
#>  80  0.2      0          5       0      
#>  120 0.4      0.4        5       0      
#>  160 0.4      0          5       0      
#> 
#> Outcome: continuous
attr(plan, "recommended_n")
#> [1] NA
```
