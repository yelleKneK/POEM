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

- outcome:

  Outcome type: `"continuous"` (linear model), `"binary"` (logistic
  model), or `"count"` (Poisson log-link model).

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

A tidy `data.frame` of class `poem_tbl`. See
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md)
for the row schema and the attributes. The identified active mediators
are in the `"active_mediators"` attribute and the print footer; the
tuning parameter HBIC chose is in the `"tuning"` attribute and the
footer.

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
identifies which individual mediators are active, with familywise error
rate control under `method = "Bonferroni"` or false discovery rate
control under `"BH"` / `"BY"`.

## When to use POEM

POEM is for the global question, “is there any active mediator among
many candidates?”, followed by the identification of the active ones. It
assumes a linear or generalized linear outcome model, a sparse set of
active mediators, standardized inputs, and complete data. Some practical
limits, seen on simulated data during the package's release audit:

- Sample size. The penalized selection needs \\n\\ large relative to
  \\\log p\\ and to the signal. With \\n\\ near 100 and \\p\\ in the
  thousands under the contrasting pattern the selection keeps a single
  mediator and the power-enhanced test has nothing to add. The default
  grid's floor rises as \\n\\ falls (see
  [`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md)),
  which protects the identified set but costs global power at small
  designs.

- Strongly correlated mediators. With an autoregressive correlation of
  0.9 among neighboring mediators the selection keeps one representative
  of a correlated block (an active mediator was recovered in 3 of 10
  seeds, against 10 of 10 at a correlation of 0.5), so the identified
  set is unstable from sample to sample even when the global test
  rejects.

- Rare binary outcomes. With about 2 to 3 percent events at \\n = 300\\
  the penalized logistic fit selects nothing at any grid value; the fit
  is empty and the function says so.

- Overdispersed counts. The count model is Poisson; under negative
  binomial overdispersion its size was not inflated at \\n = 200\\, \\p
  = 100\\, but about a tenth of the fits were empty. The count path is
  also the slowest, 15 to 30 times the continuous one.

- Missing values are not handled; supply complete cases.

For a single mediator, or a few mediators fit as one structural model
(indirect effects with confidence intervals, likelihood ratio tests of
arbitrary indirect effects, moderated mediation), the DMAR package is
the tool; POEM's contribution is the high-dimensional global test. The
other high-dimensional tests the article benchmarks against (HILMA,
GlobalTest, HDMT, DACT) live in their own packages; the vignette
`poem-vs-competitors` shows how to run them beside a POEM fit.

## References

Yu, X., and Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*. (The article these methods implement.)

Fan, J., Liao, Y., and Yao, J. (2015). Power enhancement in
high-dimensional cross-sectional tests. *Econometrica, 83*(4),
1497–1541. [doi:10.3982/ECTA12749](https://doi.org/10.3982/ECTA12749)

Guo, X., Li, R., Liu, J., and Zeng, M. (2022). High-dimensional
mediation analysis for selecting DNA methylation loci mediating
childhood trauma and cortisol stress reactivity. *Journal of the
American Statistical Association, 117*(539), 1110–1121.
[doi:10.1080/01621459.2022.2053136](https://doi.org/10.1080/01621459.2022.2053136)

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery
rate: A practical and powerful approach to multiple testing. *Journal of
the Royal Statistical Society, Series B, 57*(1), 289–300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

Benjamini, Y., and Yekutieli, D. (2001). The control of the false
discovery rate in multiple testing under dependency. *The Annals of
Statistics, 29*(4), 1165–1188.
[doi:10.1214/aos/1013699998](https://doi.org/10.1214/aos/1013699998)

## See also

The outcome-specific workers
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md);
[`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md)
for the default tuning grid;
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md)
to generate data;
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
to reproduce the article's size and power studies.

Other mediation tests:
[`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md),
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
d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
                             outcome = "continuous", c1 = 1)
pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
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
