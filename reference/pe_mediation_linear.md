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
  lambda_grid = seq(0.05, 1, length.out = 20),
  lambda_grid_reduced = lambda_grid
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
  columns; \\p\\ may exceed \\n\\.

- Z:

  Optional numeric matrix of confounders with \\n\\ rows. Default `NULL`
  (no confounders).

- method:

  Multiplicity adjustment used to screen individual mediators in the PE
  component: `"Bonferroni"` (family-wise error rate control, the
  default), `"BH"`, or `"BY"` (false discovery rate control; `"BH"`
  assumes independence or positive dependence among mediators, `"BY"` is
  valid under arbitrary dependence).

- scale:

  Logical; if `TRUE` (default) the exposures, mediators, and confounders
  are standardized (and a continuous outcome is centered) before
  fitting, as the method assumes.

- error_level:

  Target error rate for the mediator-screening step, in \\(0, 1)\\.
  Interpreted as the family-wise error rate when `method = "Bonferroni"`
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
  fit; the value minimizing the high-dimensional BIC is chosen. Default
  `seq(0.05, 1, length.out = 20)`.

- lambda_grid_reduced:

  Numeric vector of candidate tuning parameters for the reduced-model
  fit (continuous outcome only). Defaults to `lambda_grid`.

## Value

A tidy `data.frame` of class `poem_tbl` with rows `stat_hdmm` and
`pval_hdmm` (the benchmark Wald test), `stat_pe`, `j_pe`, and `pval_pe`
(the power-enhanced test and its PE component), `total_indirect_effect`
(one row per exposure when `q > 1`), `n_active_mediators`, `df`,
`n_candidate_mediators`, and `n_observations`. The identified active
mediators (their column indices in `M`) are returned in the
`"active_mediators"` attribute and shown in the print footer.

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
# near powerless here; the PE test detects the active mediators.
d <- simulate_mediation_data(n = 150, p = 80, pattern = "contrasting",
                             outcome = "continuous", c1 = 0.5)
pe_mediation_linear(d$X, d$Y, d$M)
#>  term                  value 
#>  stat_hdmm             0.1342
#>  pval_hdmm             0.7141
#>  stat_pe               0.1342
#>  j_pe                  0     
#>  pval_pe               0.7141
#>  total_indirect_effect 0.0292
#>  total_indirect_lower  -0.127
#>  total_indirect_upper  0.1854
#>  n_active_mediators    0     
#>  df                    1     
#>  n_candidate_mediators 80    
#>  n_observations        150   
#> 
#> Outcome model: continuous (linear)
#> Active mediators identified: none
```
