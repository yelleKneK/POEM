# Monte Carlo study of individual-mediator identification (FWER / FDR)

Evaluates how well the power-enhanced screen recovers the *individual*
active mediators, the experiment reported in the manuscript's
supplement. For each signal-strength scale \\c_1\\ on a grid it
simulates many data sets, identifies the active set under each
multiplicity method, and scores those selections against the known
truth, returning the empirical family-wise error rate, false discovery
rate, precision, and recall. This complements
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md),
which evaluates the *global* test rather than which mediators are
flagged.

## Usage

``` r
pe_identification_study(
  n,
  p,
  outcome = c("continuous", "binary", "count"),
  pattern = c("homogeneous", "contrasting", "heterogeneous"),
  c1_grid = seq(0, 1, by = 0.25),
  c2 = 0.5,
  n_rep = 200,
  methods = c("Bonferroni", "BH", "BY"),
  error_level = 0.05,
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
  `seq(0, 1, by = 0.25)`. Include 0 to estimate the null family-wise
  error rate.

- c2:

  Direct effect used in the simulation. Default 0.5.

- n_rep:

  Number of Monte Carlo replications per grid point. Default 200. The
  manuscript uses 1000.

- methods:

  Multiplicity methods to evaluate, any of `"Bonferroni"` (family-wise
  error rate), `"BH"`, and `"BY"` (false discovery rate). Default all
  three.

- error_level:

  Target error rate for the mediator-screening step (the FWER level for
  Bonferroni, the FDR level for BH and BY). Default 0.05.

- outcome_args:

  A list of further arguments forwarded to
  [`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md)
  (for example `rho`, `q`, `d`, `tau`, `alpha_m`). Default empty.

- cores:

  Number of CPU cores for the replications. Values above 1 fork via the
  base parallel package (Unix only); a seeded parallel run is
  reproducible across runs at the same `cores` but need not match a
  serial run. Default 1.

- seed:

  Optional integer seed, set locally with the caller's random number
  generator state restored on exit.

- progress:

  Logical; if `TRUE`, print a line per grid point. Default `FALSE`.

## Value

A tidy `data.frame` of class `poem_tbl` with one row per (\\c_1\\,
method) combination and columns `c1`, `method`, `fwer` (empirical
family-wise error rate, the proportion of replications with at least one
false positive), `fdr` (mean false discovery proportion), `precision`,
`recall`, and `n_valid` (replications that fit successfully). The
simulation settings are recorded in attributes.

## Details

At \\c_1 = 0\\ the exposure-on-mediator paths are all zero, so there are
no truly active mediators and any selection is a false positive; the
family-wise error rate column is then the probability of selecting any
mediator (the empirical FWER under the null). Recall is undefined when
there are no active mediators (at \\c_1 = 0\\) and is reported as `NA`
there; precision is `NA` for a replication that selects nothing and is
averaged over the replications that select at least one mediator.

## References

Yu, X., and Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

## See also

[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
for the global test,
[`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md)
for the per-method active set of a single fit.

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEM/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEM/reference/mediation_ar1_cov.md),
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
# \donttest{
# A small, fast study. Raise n, p, and n_rep toward the manuscript's
# settings for a publication-scale study.
set.seed(113)
pe_identification_study(n = 150, p = 60, outcome = "continuous",
                        pattern = "contrasting", c1_grid = c(0, 0.5, 1),
                        n_rep = 20)
#>  c1  method     fwer fdr     precision recall n_valid
#>  0   Bonferroni 0.05 0.05    0         <NA>   20     
#>  0   BH         0.05 0.05    0         <NA>   20     
#>  0   BY         0.05 0.05    0         <NA>   20     
#>  0.5 Bonferroni 0    0       1         0.1125 20     
#>  0.5 BH         0    0       1         0.15   20     
#>  0.5 BY         0    0       1         0.1    20     
#>  1   Bonferroni 0    0       1         0.5125 20     
#>  1   BH         0.15 0.03917 0.9565    0.5625 20     
#>  1   BY         0    0       1         0.4875 20     
#> 
#> Outcome model: continuous
# }
```
