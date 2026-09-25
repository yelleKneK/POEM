# Power-enhanced test for high-dimensional mediation

Tests whether any mediator is active among a large, possibly
intercorrelated set of candidate mediators, for a continuous, binary, or
count outcome. This is the main entry point of POEM: name the outcome
type and it dispatches to the appropriate model, returning one tidy
table that reports the benchmark Wald test on the total indirect effect
side by side with the power-enhanced (PE) test, plus the set of
individual mediators identified as active.

## Usage

``` r
pe_mediation(
  X,
  Y,
  M,
  Z = NULL,
  outcome = c("continuous", "binary", "count"),
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

- outcome:

  Outcome type: `"continuous"` (linear model), `"binary"` (logistic
  model), or `"count"` (Poisson log-link model).

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

A tidy `data.frame` of class `poem_tbl`. See
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md)
for the row schema. The identified active mediators are in the
`"active_mediators"` attribute and the print footer.

## Details

Standard high-dimensional mediation tests target the *total* indirect
effect \\\beta = \Gamma_x \alpha_m\\, the sum of every mediator's
individual indirect effect. When some indirect effects are positive and
others negative they can cancel, making \\\beta = 0\\ even though active
mediators exist, and a test built on \\\beta\\ is then powerless. The
power-enhanced test adds a component \$\$J_m = \sqrt{p} \sum\_{i=1}^{q}
\sum\_{j \in \hat{S}} \left\|
\frac{\hat\alpha\_{m,j}}{\hat\sigma\_{m,j}} \right\| \left\|
\frac{\hat\Gamma\_{x,i,j}}{\hat\sigma\_{\Gamma,i,j}} \right\|
\mathbf{1}\\\left\\ \max(p\_{1,j}, p\_{2,i,j}) \<
\frac{\alpha\_{\mathrm{lvl}}}{s q \log\log n} \right\\,\$\$ a sum over
the selected mediators \\\hat S\\ of the product of the (standardized)
mediator-on-outcome and exposure-on-mediator statistics, kept only for
pairs where both paths are individually significant. Because \\J_m\\
accumulates magnitudes, opposite-signed indirect effects reinforce
rather than cancel, so the test \\M\_{PE} = S_n + J_m\\ stays powerful
under heterogeneous and contrasting mediation while keeping the
chi-square reference distribution (and hence the Type I error rate) of
the benchmark test \\S_n\\ under the global null. The indicator also
identifies which individual mediators are active, with family-wise error
rate control under `method = "Bonferroni"` or false discovery rate
control under `"BH"` / `"BY"`.

## References

Yu, X., and Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*. (The article these methods implement.)

Fan, J., Liao, Y., and Yao, J. (2015). Power enhancement in
high-dimensional cross-sectional tests. *Econometrica, 83*(4),
1497–1541.

Guo, X., Li, R., Liu, J., and Zeng, M. (2022). High-dimensional
mediation analysis for selecting DNA methylation loci mediating
childhood trauma and cortisol stress reactivity. *Journal of the
American Statistical Association, 117*(539), 1110–1121.

## See also

The outcome-specific workers
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md);
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md)
to generate data;
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
to reproduce the manuscript's size and power studies.

Other mediation tests:
[`pe_mediate()`](https://yelleknek.github.io/POEM/reference/pe_mediate.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
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
# Contrasting mediation: active mediators whose effects cancel, so the
# total indirect effect is zero. Compare the two p-values: the benchmark
# Wald test (pval_hdmm) is large while the PE test (pval_pe) is tiny.
d <- simulate_mediation_data(n = 150, p = 80, pattern = "contrasting",
                             outcome = "continuous", c1 = 0.5)
pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
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
