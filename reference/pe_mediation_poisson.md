# Power-enhanced mediation test for a count outcome

Tests the global null hypothesis of no active mediator for a count
outcome modeled with Poisson (log-link) regression, reporting both the
benchmark Wald test on the total indirect effect and the power-enhanced
(PE) test. This is the count-outcome worker behind
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md).

## Usage

``` r
pe_mediation_poisson(
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

A tidy `data.frame` of class `poem_tbl` with the same schema as
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md).

## Details

The outcome follows a Poisson mediation model, \\Y \mid M, X, Z \sim
\mathrm{Poisson}(\exp(\alpha_m' M + \alpha_x' X + \alpha_z' Z))\\, with
a linear mediator model as in
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md).
Estimation, tuning, and the power enhancement construction follow the
generalized-outcome framework of Guo et al. (2024); the penalized
Poisson fit uses a glmnet lasso warm start followed by a weighted local
linear approximation step. The total indirect effect \\\beta = \Gamma_x
\alpha_m\\ is on the log-mean (link) scale, not the count scale.

## References

Guo, X., Li, R., Liu, J., and Zeng, M. (2024). Estimations and tests for
generalized mediation models with high-dimensional potential mediators.
*Journal of Business & Economic Statistics, 42*(1), 243–256.
[doi:10.1080/07350015.2023.2174548](https://doi.org/10.1080/07350015.2023.2174548)

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md).

Other mediation tests:
[`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md),
[`pe_mediate()`](https://yelleknek.github.io/POEM/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md),
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md),
[`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md),
[`summary.poem_tbl()`](https://yelleknek.github.io/POEM/reference/summary.poem_tbl.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
                             outcome = "count", c1 = 0.5, c2 = 0.4)
pe_mediation_poisson(d$X, d$Y, d$M)
#>  term                  value     
#>  stat_hdmm             0.0002879 
#>  pval_hdmm             0.9865    
#>  stat_pe               545.4     
#>  j_pe                  545.4     
#>  pval_pe               < 0.0001  
#>  total_indirect_effect -0.0007606
#>  total_indirect_lower  -0.08861  
#>  total_indirect_upper  0.08709   
#>  n_active_mediators    2         
#>  df                    1         
#>  n_candidate_mediators 60        
#>  n_observations        200       
#> 
#> Outcome model: count (Poisson)
#> Active mediators identified (2): 4, 5
#> Tuning parameter (HBIC): lambda = 0.2 from 20 values in [0.05, 1]
```
