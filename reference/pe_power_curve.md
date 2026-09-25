# Monte Carlo size and power curve for the PE mediation tests

Reproduces the article's simulation studies: for each value of the
signal-strength scale \\c_1\\ on a grid, it simulates many data sets,
applies both the benchmark Wald test and the power-enhanced test, and
returns their empirical rejection rates. At \\c_1 = 0\\ (the global null
of no mediation) the rejection rate estimates the Type I error rate; at
nonzero \\c_1\\ it estimates power. The headline finding is visible
directly in the table: under a contrasting pattern the benchmark test
stays near its nominal size as \\\|c_1\|\\ grows while the
power-enhanced test climbs toward one.

## Usage

``` r
pe_power_curve(
  n,
  p,
  outcome = c("continuous", "binary", "count"),
  pattern = c("homogeneous", "contrasting", "heterogeneous"),
  c1_grid = seq(0, 1, by = 0.25),
  c2 = 0.5,
  n_rep = 100,
  alpha_level = 0.05,
  method = c("Bonferroni", "BH", "BY"),
  error_level = 0.05,
  lambda_grid = NULL,
  lambda_grid_reduced = NULL,
  outcome_args = list(),
  cores = 1L,
  seed = NULL,
  progress = FALSE
)
```

## Arguments

- n, p:

  Number of observations and candidate mediators per data set.

- outcome:

  Outcome type passed to
  [`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md):
  `"continuous"`, `"binary"`, or `"count"`.

- pattern:

  Mediation pattern: `"homogeneous"` or `"contrasting"`
  (`"heterogeneous"` is a synonym for `"contrasting"`).

- c1_grid:

  Numeric vector of signal-strength scales to sweep. Default
  `seq(0, 1, by = 0.25)`. Include 0 to estimate the Type I error rate;
  the article also uses negative values.

- c2:

  Direct effect used in the simulation. Default 0.5.

- n_rep:

  Number of Monte Carlo replications per grid point. Default 100. The
  article uses 1000.

- alpha_level:

  Significance level for the global test. A replication counts as a
  rejection when its p-value is at most `alpha_level`. Default 0.05.

- method, error_level:

  Passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md):
  the multiplicity method and target error rate for the
  mediator-screening step.

- lambda_grid, lambda_grid_reduced:

  Passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md):
  the tuning grids for the penalized fit. Default `NULL`, the family's
  default grid (see
  [`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md)).

- outcome_args:

  A list of further arguments forwarded to
  [`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md)
  (for example `rho`, `q`, `d`, `tau`, `alpha_m`). Default empty. The
  design arguments this function sets itself (`n`, `p`, `outcome`,
  `pattern`, `c1`, `c2`) and `seed` may not appear here; a `seed` inside
  `outcome_args` would make every replication draw the same data.

- cores:

  Number of CPU cores for the Monte Carlo replications. Values above 1
  fork via the base parallel package (Unix only). With parallel cores,
  reproducibility uses parallel-safe streams, so a seeded parallel run
  is reproducible across runs at the same `cores` but need not match a
  serial run. Default 1.

- seed:

  Optional integer seed, set locally with the caller's random number
  generator state restored on exit.

- progress:

  Logical; if `TRUE`, print a line per grid point as it completes.
  Default `FALSE`.

## Value

A tidy `data.frame` of class `poem_tbl` with one row per grid point and
columns `c1`, `rejection_hdmm`, and `rejection_pe` (the empirical
rejection rates of the benchmark and power-enhanced tests), `n_valid`
(replications whose fit succeeded, the denominator of the rates), and
`n_empty` (replications whose penalized fit selected no mediator,
counted as non-rejections). The simulation settings are recorded in the
attributes `outcome`, `pattern`, `n`, `p`, `n_rep`, `alpha_level`,
`error_level`, and `lambda_grid`.

## Details

A rejection rate from `n_rep` replications carries a Monte Carlo
standard error of about \\\sqrt{r(1 - r) / n\_{rep}}\\; at 100
replications a rate near 0.05 is known to about 0.02 and a rate near 0.5
to about 0.05. The article uses 1000 replications. A replication whose
penalized fit selects no mediator counts as a non-rejection (both
p-values are 1), the article's convention; the number of such
replications is reported in `n_empty`. A replication whose fit fails
outright is dropped from the denominator and reported in `n_valid`, with
a warning that quotes the first error.

The tuning grid matters for what these curves show. The article's
figures were produced with per-setting grids (see
[`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md));
the default here is that grid rescaled to `n` and `p`, so the article's
settings reproduce its curves within Monte Carlo error. Pass
`lambda_grid` to study another grid.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md),
[`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md).

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEM/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEM/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEM/reference/pe_identification_study.md),
[`pe_simulation_study()`](https://yelleknek.github.io/POEM/reference/pe_simulation_study.md),
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
pe_power_curve(n = 200, p = 60, outcome = "continuous",
               pattern = "contrasting", c1_grid = c(0, 0.5, 1),
               n_rep = 5)
#>  c1  rejection_hdmm rejection_pe n_valid n_empty
#>  0   0              0            5       0      
#>  0.5 0              0.2          5       0      
#>  1   0.2            0.4          5       0      
#> 
#> Outcome: continuous
```
