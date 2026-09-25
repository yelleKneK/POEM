# Summarize a POEM table

For a grouped table that compares the benchmark and power-enhanced tests
(any `poem_tbl` with a `group` column, a `pval_hdmm` column, and a
power-enhanced p-value column (`pval_pe` or `pval_pe_bonferroni`), such
as the output of
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md)),
[`summary()`](https://rdrr.io/r/base/summary.html) digests it across
groups: how many groups each test detects mediation in, and in
particular the groups where the power-enhanced test detects mediation
that the benchmark misses, together with the most frequently flagged
mediators. For any other `poem_tbl` it falls back to the ordinary
data-frame summary.

## Usage

``` r
# S3 method for class 'poem_tbl'
summary(object, alpha_level = 0.05, ...)
```

## Arguments

- object:

  A `poem_tbl`.

- alpha_level:

  Significance level for counting a group as a detection. Default 0.05.

- ...:

  Ignored.

## Value

For a grouped comparison table, an object of class `summary.poem_tbl` (a
list, printed by its own method) with elements `outcome`, `alpha_level`,
`n_groups`, `n_with_data`, `n_hdmm`, `n_pe` (groups detected by each
test), `pe_only` (a `data.frame` of the groups the PE test detects but
the benchmark does not), and `top_mediators` (a frequency table of
flagged mediators). For any other `poem_tbl`, the value of
[`summary.data.frame()`](https://rdrr.io/r/base/summary.html).

## See also

[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md),
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md).

Other mediation tests:
[`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md),
[`pe_mediate()`](https://yelleknek.github.io/POEM/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md),
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md),
[`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# The global and income-group infant-mortality models on a 25-point
# tuning grid, which keeps the example fast; the article's tables use
# the default 100-point grid.
res <- WHO_mediation_analysis("imr", groupings = c("global", "income"),
                              lambda_grid = seq(0.1, 10, length.out = 25))
summary(res)
#> POEM comparison: IMR
#>   groups: 5 (5 with data), alpha_level = 0.05
#>   detected by benchmark (HDMM): 2   by power-enhanced (PE): 3
#>   PE detects mediation in 1 group the benchmark misses:
#>  group n_countries pval_hdmm  pval_pe active_mediators
#>    Low          16     0.432 7.82e-49     pvtd_usd2021
#>   most-flagged mediators:
#>     gge_gdp            1
#>     oops_che           1
#>     pvtd_usd2021       1
#>     shi_che            1
```
