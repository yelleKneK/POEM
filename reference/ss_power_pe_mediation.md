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
  alpha = 0.05,
  method = c("Bonferroni", "BH", "BY"),
  error_level = 0.05,
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

- alpha:

  Significance level. Default 0.05.

- method, error_level:

  Passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md).

- cores:

  Number of CPU cores (Unix forking). Default 1.

- seed:

  Optional integer seed, set locally and restored on exit.

## Value

A tidy `data.frame` (class `poem_tbl`) with one row per candidate sample
size and columns `n`, `power_pe`, `power_hdmm`, and `n_valid`. The
smallest `n` reaching `target_power` for the PE test is stored in the
`"recommended_n"` attribute (`NA` if no grid value reaches it).

## Details

This is a simulation-based procedure, not a closed-form power formula:
no analytic power expression exists for the PE test, so sample size is
instead determined by direct Monte Carlo simulation. Power is estimated
by simulating `n_rep` data sets at each candidate `n` with
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md)
and recording how often the test rejects at level `alpha`. Because the
estimate is a Monte Carlo proportion, it carries simulation error of
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
# \donttest{
set.seed(113)
plan <- ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
                              p = 50, n_grid = c(80, 120, 160), n_rep = 30)
plan
#>  n   power_pe power_hdmm n_valid
#>  80  0.5333   0.1333     30     
#>  120 0.8667   0.06667    30     
#>  160 0.9667   0.1        30     
#> 
#> Outcome model: continuous
attr(plan, "recommended_n")
#> [1] 120
# }
```
