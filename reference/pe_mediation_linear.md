# Power-enhanced mediation test for a continuous outcome

Tests the global null hypothesis of no active mediator among a
high-dimensional set of candidate mediators, for a continuous
(linear-model) outcome, and reports both the benchmark Wald test on the
total indirect effect and the power-enhanced (PE) test that remains
powerful when individual mediation effects are heterogeneous or
contrasting. This is the continuous-outcome worker behind
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md);
call it directly when the outcome is continuous.

## Usage

``` r
pe_mediation_linear(
  X,
  Y,
  M,
  Z = NULL,
  method = c("Bonferroni", "BH", "BY"),
  scale = TRUE,
  error_level = 0.05,
  conf_level = 0.95,
  report_all_methods = FALSE,
  drop_constant = FALSE,
  lambda_grid = NULL,
  lambda_grid_reduced = NULL
)
```

## Arguments

- X:

  Numeric matrix of exposures with \\n\\ rows and \\q\\ columns (\\q\\
  is usually 1). A numeric vector is treated as a single exposure.

- Y:

  Numeric outcome vector of length \\n\\. Continuous for
  `outcome = "continuous"`, coded 0/1 for `"binary"`, and nonnegative
  integer counts for `"count"`.

- M:

  Numeric matrix of candidate mediators with \\n\\ rows and \\p\\
  columns; \\p\\ may exceed \\n\\. Column names, when present, label the
  active mediators in the print footer and the mediator table; without
  them mediators are reported by column position.

- Z:

  Optional numeric matrix of confounders with \\n\\ rows. Default `NULL`
  (no confounders).

- method:

  Multiplicity adjustment used to screen individual mediators in the PE
  component: `"Bonferroni"` (familywise error rate control, the
  default), `"BH"`, or `"BY"` (false discovery rate control; `"BH"`
  assumes independence or positive dependence among mediators, `"BY"` is
  valid under arbitrary dependence).

- scale:

  Logical; if `TRUE` (default) the exposures, mediators, and confounders
  are standardized (and a continuous outcome is centered) before
  fitting, as the method assumes.

- error_level:

  Target error rate for the mediator-screening step, in \\(0, 1)\\.
  Interpreted as the familywise error rate when `method = "Bonferroni"`
  and as the false discovery rate when `method = "BH"` or `"BY"`.
  Default 0.05. This is distinct from the significance level used to
  test the global null (which is the user's choice when reading
  `pval_pe`); the test is insensitive to `error_level` and remains valid
  for any fixed value in \\(0, 1)\\.

- conf_level:

  Confidence level for the interval on the total indirect effect.
  Default 0.95.

- report_all_methods:

  Logical; if `TRUE`, the active-mediator set is computed for all three
  multiplicity methods (Bonferroni, BH, BY) off the same fit and
  recorded for comparison, retrievable with
  [`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md).
  The `method` argument still drives the primary result. Default
  `FALSE`.

- drop_constant:

  Logical; how to handle constant (zero-variance) mediator columns,
  which cannot be standardized or carry signal. If `FALSE` (the default)
  they are an error; if `TRUE` they are dropped with a warning and the
  remaining mediators are renumbered back to their original column
  positions in the reported active set. Useful for small subgroups where
  some mediators happen to be constant.

- lambda_grid:

  Numeric vector of candidate SCAD tuning parameters for the penalized
  mediator fit; the fit is repeated at each value and the one minimizing
  the high-dimensional BIC (HBIC) is kept. Default `NULL` uses, for a
  continuous outcome,
  [`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md):
  the grid the article's simulations used at their design, rescaled to
  this `n` and `p` by the rate \\\sqrt{\log p / n}\\ the theory requires
  of the tuning parameter; for a binary or count outcome it uses the
  review-era grid `seq(0.05, 1, length.out = 20)`. The grid's lower end
  matters most, because the HBIC minimum often sits there and because it
  decides whether the identified set keeps its error-rate guarantee; see
  [`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md)
  for the trade-off with power, measured. The value chosen is reported
  in the `"tuning"` attribute and the print footer.

- lambda_grid_reduced:

  Numeric vector of candidate tuning parameters for the reduced-model
  fit that the continuous-outcome Wald test refits (unused for binary
  and count outcomes). Default `NULL` uses
  [`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md)
  with `model = "reduced"` when `lambda_grid` is also `NULL` (the
  article's reduced-model grid, rescaled), and otherwise `lambda_grid`
  itself.

## Value

A tidy `data.frame` of class `poem_tbl` with rows `stat_hdmm` and
`pval_hdmm` (the benchmark Wald test), `stat_pe`, `j_pe`, and `pval_pe`
(the power-enhanced test and its PE component), `total_indirect_effect`
with `total_indirect_lower` and `total_indirect_upper` (the estimate and
a Wald interval at `conf_level`, the estimate plus or minus a standard
normal quantile times its standard error; one row of each per exposure
when `q > 1`, suffixed with the exposure's column name when `X` has
column names and with `_1`, `_2`, ... otherwise), `n_active_mediators`,
`df` (the chi-square degrees of freedom, equal to `q`),
`n_candidate_mediators`, and `n_observations`. Attributes:
`"active_mediators"` (the column positions in `M` identified as active;
the print footer shows their column names when `M` has them),
`"mediator_names"` and `"exposure_names"` (the column names of `M` and
`X`, or `NULL` when unnamed), `"mediator_table"` (the per-mediator
screening statistics, see
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md)),
`"method"`, `"error_level"`, `"conf_level"`, `"outcome"`, `"tuning"` (a
list with the grid searched, the HBIC-selected `lambda_selected`, and
`at_lower_end`, `TRUE` when the selection sat at the grid's smallest
value; for a continuous outcome also the reduced model's grid and
selection), `"empty_fit"`, and, with `report_all_methods = TRUE`,
`"selection_by_method"` and `"pe_by_method"` (see
[`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md)).

When the penalized fit selects no mediator at any grid value there is
nothing to test: the function warns (a condition of class
`poem_empty_fit`), sets `"empty_fit"` to `TRUE`, and returns the same
rows with both statistics 0, both p-values 1 (the article's Monte Carlo
convention, under which an empty selection is a non-rejection), a total
indirect effect of 0, and `NA` interval limits. A p-value of 1 from such
a fit is not evidence for the null; a grid reaching smaller values may
select mediators.

## Details

The model is the linear mediation pair \$\$Y = \alpha_m' M + \alpha_x'
X + \alpha_z' Z + \varepsilon_y, \qquad M = \Gamma_x' X + \Gamma_z' Z +
\varepsilon_m,\$\$ with total indirect effect \\\beta = \Gamma_x
\alpha_m\\. The mediator coefficients \\\alpha_m\\ are estimated by
partial penalized least squares with a SCAD penalty, the tuning
parameter chosen by the high-dimensional BIC. The benchmark statistic is
the Wald statistic \\S_n = n \hat\beta' \hat\Sigma\_\beta^{-1}
\hat\beta\\ of Guo et al. (2022); the power-enhanced statistic adds the
component \\J_m\\ that accumulates the marginal signal from each
selected mediator (see
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
for the rationale and the formula). Both are referred to a chi-square
distribution with \\q\\ (the number of exposures) degrees of freedom
under the global null.

## References

Guo, X., Li, R., Liu, J., and Zeng, M. (2022). High-dimensional
mediation analysis for selecting DNA methylation loci mediating
childhood trauma and cortisol stress reactivity. *Journal of the
American Statistical Association, 117*(539), 1110–1121.
[doi:10.1080/01621459.2022.2053136](https://doi.org/10.1080/01621459.2022.2053136)

Fan, J., and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *Journal of the American
Statistical Association, 96*(456), 1348–1360.
[doi:10.1198/016214501753382273](https://doi.org/10.1198/016214501753382273)

Wang, L., Kim, Y., and Li, R. (2013). Calibrating nonconvex penalized
regression in ultra-high dimension. *The Annals of Statistics, 41*(5),
2505–2536. [doi:10.1214/13-AOS1159](https://doi.org/10.1214/13-AOS1159)

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
for the front end,
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md)
and
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md)
for binary and count outcomes, and
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md)
to generate data for trying the method.

Other mediation tests:
[`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md),
[`pe_mediate()`](https://yelleknek.github.io/POEM/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md),
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md),
[`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md),
[`summary.poem_tbl()`](https://yelleknek.github.io/POEM/reference/summary.poem_tbl.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
# A contrasting setting: two active mediators whose indirect effects
# cancel, so the total indirect effect is zero. The Wald test is
# near powerless here; the PE test rejects and identifies active
# mediators.
d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
                             outcome = "continuous", c1 = 1)
pe_mediation_linear(d$X, d$Y, d$M)
#>  term                  value   
#>  stat_hdmm             0.03881 
#>  pval_hdmm             0.8438  
#>  stat_pe               174.2   
#>  j_pe                  174.1   
#>  pval_pe               < 0.0001
#>  total_indirect_effect 0.01054 
#>  total_indirect_lower  -0.0943 
#>  total_indirect_upper  0.1154  
#>  n_active_mediators    1       
#>  df                    1       
#>  n_candidate_mediators 60      
#>  n_observations        200     
#> 
#> Outcome model: continuous (linear)
#> Active mediators identified (1): 2
#> Tuning parameter (HBIC): lambda = 0.199 from 20 values in [0.199, 0.388] (the grid's lower end)
```
