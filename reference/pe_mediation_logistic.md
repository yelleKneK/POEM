# Power-enhanced mediation test for a binary outcome

Tests the global null hypothesis of no active mediator for a binary
(0/1) outcome modeled with logistic regression, reporting both the
benchmark Wald test on the total indirect effect and the power-enhanced
(PE) test. This is the binary-outcome worker behind
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md).

## Usage

``` r
pe_mediation_logistic(
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

A tidy `data.frame` of class `poem_tbl` with the same schema as
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md).

## Details

The outcome follows a logistic mediation model, \\\mathrm{logit}\\P(Y =
1 \mid M, X, Z) = \alpha_m' M + \alpha_x' X + \alpha_z' Z\\, with a
linear mediator model as in
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md).
The mediator coefficients are estimated by partial penalized likelihood
with a SCAD penalty (a coordinate-descent fit that leaves the exposures
and confounders unpenalized), the tuning parameter chosen by the
high-dimensional BIC. On the logit (link) scale the total indirect
effect \\\beta = \Gamma_x \alpha_m\\ and the power enhancement
construction are exactly as in the continuous case (Guo et al., 2024).
The total indirect effect here is on the log-odds scale, not the
probability scale.

## References

Guo, X., Li, R., Liu, J., Zeng, M., and Wang, H. (2024). Statistical
inference for high-dimensional generalized mediation analysis. (Working
paper.)

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md).

Other mediation tests:
[`pe_mediate()`](https://yelleknek.github.io/POEM/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md),
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md),
[`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md),
[`summary.poem_tbl()`](https://yelleknek.github.io/POEM/reference/summary.poem_tbl.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
                             outcome = "binary", c1 = 0.5, c2 = 1)
pe_mediation_logistic(d$X, d$Y, d$M)
#>  term                  value  
#>  stat_hdmm             1.842  
#>  pval_hdmm             0.1748 
#>  stat_pe               1.842  
#>  j_pe                  0      
#>  pval_pe               0.1748 
#>  total_indirect_effect -0.282 
#>  total_indirect_lower  -0.6894
#>  total_indirect_upper  0.1253 
#>  n_active_mediators    0      
#>  df                    1      
#>  n_candidate_mediators 60     
#>  n_observations        200    
#> 
#> Outcome model: binary (logistic)
#> Active mediators identified: none
```
