# Active mediators identified under each multiplicity method

Summarizes which mediators a
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
fit identifies as active. When the fit was produced with
`report_all_methods = TRUE`, this reports the active set under each of
the three multiplicity methods (Bonferroni for family-wise error rate
control, Benjamini-Hochberg and Benjamini-Yekutieli for false discovery
rate control), so their conservativeness can be compared on the same
fit. Otherwise it reports the single method that was used.

## Usage

``` r
pe_selection(fit)
```

## Arguments

- fit:

  A result from
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
  or one of its workers.

## Value

A `data.frame` with one row per multiplicity method and columns
`method`, `n_active` (number of active mediators), `j_pe` (the power
enhancement component \\J_m\\ under that method's screen), `stat_pe` and
`pval_pe` (the resulting power-enhanced statistic \\M\_{PE}\\ and its
p-value), and `active_mediators` (their column indices in `M`,
comma-separated, or `"none"`). Because each multiplicity method screens
a different active set into \\J_m\\, the `pval_pe` column gives the PE /
PE_BH / PE_BY global p-values side by side, as the manuscript's extended
data-analysis tables report them.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
(its `report_all_methods` argument),
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md).

Other mediation tests:
[`pe_mediate()`](https://yelleknek.github.io/POEM/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md),
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md),
[`summary.poem_tbl()`](https://yelleknek.github.io/POEM/reference/summary.poem_tbl.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                             pattern = "contrasting", c1 = 1)
fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                    report_all_methods = TRUE)
pe_selection(fit)
#>       method n_active      j_pe   stat_pe       pval_pe active_mediators
#> 1 Bonferroni        3  864.1774  864.4981 5.124789e-190          2, 3, 4
#> 2         BH        4 1105.1162 1105.4370 2.173658e-242      2, 3, 4, 33
#> 3         BY        3  864.1774  864.4981 5.124789e-190          2, 3, 4
```
