# Monte Carlo study of individual-mediator identification (FWER / FDR)

Evaluates how well the power-enhanced screen recovers the *individual*
active mediators, the experiment reported in the article's supplement.
For each signal-strength scale \\c_1\\ on a grid it simulates many data
sets, identifies the active set under each multiplicity method, and
scores those selections against the known truth, returning the empirical
familywise error rate, false discovery rate, precision, and recall. This
complements
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
  truth = c("outcome_effect", "mediation"),
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
  `seq(0, 1, by = 0.25)`. Include 0 to estimate the null familywise
  error rate.

- c2:

  Direct effect used in the simulation. Default 0.5.

- n_rep:

  Number of Monte Carlo replications per grid point. Default 200. The
  article uses 1000.

- methods:

  Multiplicity methods to evaluate, any of `"Bonferroni"` (familywise
  error rate), `"BH"`, and `"BY"` (false discovery rate). Default all
  three.

- error_level:

  Target error rate for the mediator-screening step (the FWER level for
  Bonferroni, the FDR level for BH and BY). Default 0.05.

- truth:

  Which mediators count as truly active: `"outcome_effect"` (nonzero
  outcome coefficient, the article's convention) or `"mediation"` (both
  paths nonzero). See Details.

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
  `pattern`, `c1`, `c2`) and `seed` may not appear here.

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
method) combination and columns `c1`, `method`, `fwer`, `fdr`,
`precision`, `recall` (as defined in Details), `n_valid` (replications
whose fit succeeded), and `n_empty` (replications whose penalized fit
selected no mediator). The settings are recorded in the attributes
`outcome`, `pattern`, `n`, `p`, `n_rep`, `error_level`, `truth`, and
`lambda_grid`.

## Details

The four metrics follow the article's supplement. In each replication
the selected set is compared with the truth set; the familywise error
rate is the proportion of replications selecting at least one mediator
outside the truth set, the false discovery rate is the mean false
discovery proportion (zero when nothing is selected), precision is the
mean proportion of selected mediators that are true (scored zero when
nothing is selected), and recall is the mean proportion of true
mediators selected. Each is a Monte Carlo proportion with standard error
about \\\sqrt{r(1 - r) / n\_{rep}}\\; the article uses 1000
replications.

What counts as a true mediator is set by `truth`. The article scores a
mediator as active when its outcome coefficient \\\alpha\_{m,j}\\ is
nonzero (`truth = "outcome_effect"`, the default), which is what its
supplement tables use and what makes precision and recall defined at
\\c_1 = 0\\, where the exposure-on-mediator paths are all zero; there
the familywise error rate is the probability of selecting a mediator
with no outcome effect. Under `truth = "mediation"` a mediator is active
only when both paths are nonzero (the simulator's `active_mediators`),
so at \\c_1 = 0\\ no mediator is active, any selection is a false
positive, and recall is `NA`. For nonzero \\c_1\\ the two definitions
agree under the article's designs, whose exposure-on-mediator loadings
are all nonzero.

The tuning grid governs these rates as much as the screen does (see
[`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md));
the default reproduces the article's settings.

## References

Yu, X., and Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

## See also

[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
for the global test,
[`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md)
for the per-method active set of a single fit,
[`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md).

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
# n_rep = 5 keeps this example fast; a rate from five replications has a
# Monte Carlo standard error of up to 0.22, so read the shape, not the
# numbers. A reported study uses the article's design (n = 300, p = 500)
# and 1000 replications.
set.seed(113)
pe_identification_study(n = 200, p = 60, outcome = "continuous",
                        pattern = "contrasting", c1_grid = c(0, 1),
                        n_rep = 5)
#>  c1 method     fwer fdr precision recall n_valid n_empty
#>  0  Bonferroni 0    0   0         0      5       0      
#>  0  BH         0    0   0         0      5       0      
#>  0  BY         0    0   0         0      5       0      
#>  1  Bonferroni 0    0   0.8       0.2    5       0      
#>  1  BH         0    0   0.8       0.3    5       0      
#>  1  BY         0    0   0.8       0.25   5       0      
#> 
#> Outcome: continuous
```
