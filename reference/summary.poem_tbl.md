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
summary(object, alpha = 0.05, ...)
```

## Arguments

- object:

  A `poem_tbl`.

- alpha:

  Significance level for counting a group as a detection. Default 0.05.

- ...:

  Ignored.

## Value

For a grouped comparison table, an object of class `summary.poem_tbl` (a
list, printed by its own method) with elements `outcome`, `alpha`,
`n_groups`, `n_with_data`, `n_hdmm`, `n_pe` (groups detected by each
test), `pe_only` (a `data.frame` of the groups the PE test detects but
the benchmark does not), and `top_mediators` (a frequency table of
flagged mediators). For any other `poem_tbl`, the value of
[`summary.data.frame()`](https://rdrr.io/r/base/summary.html).

## See also

[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md),
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md).

Other mediation tests:
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
# \donttest{
res <- WHO_mediation_analysis("imr", groupings = c("global", "region"))
summary(res)
#> POEM comparison: IMR
#>   groups: 7 (7 with data), alpha = 0.05
#>   detected by benchmark (HDMM): 3   by power-enhanced (PE): 5
#>   PE detects mediation in 2 groups the benchmark misses:
#>  group n_countries pval_hdmm  pval_pe  active_mediators
#>   SEAR           5     0.309 4.13e-41    ext_usd2021_pc
#>    WPR           9     0.872 3.77e-83 chi_che, pvtd_gdp
#>   most-flagged mediators:
#>     gge_gdp            2
#>     chi_che            1
#>     ext_usd2021_pc     1
#>     pvtd_gdp           1
# }
```
