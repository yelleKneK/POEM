# Reproduce the article's empirical mediation analysis

Runs the article's real-data analysis for one health outcome across the
groupings it reports: a global model using all WHO members, then
separate models within each WHO region, within each World Bank income
group, and (optionally) within each region-by-income cell. For every
grouping it reports the benchmark Wald test and the power-enhanced test
on the total indirect effect of health expenditure, together with the
individual mediators the PE test identifies as active. This reproduces
the structure of the article's data-analysis tables (for example its
Table for infant mortality) from the shipped
[WHO_health_mediation](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md)
data.

## Usage

``` r
WHO_mediation_analysis(
  outcome = c("imr", "u5mr", "leb", "lbw", "pou"),
  groupings = c("global", "region", "income"),
  error_level = 0.05,
  lambda_grid = seq(0.1, 10, length.out = 100),
  lambda_grid_global = seq(0.1, 5, length.out = 100),
  data = WHO_health_mediation,
  cores = 1L,
  verbose = FALSE,
  full_table = FALSE
)
```

## Arguments

- outcome:

  Which health outcome to analyze: `"imr"`, `"u5mr"`, `"leb"`, `"lbw"`,
  or `"pou"`. See
  [WHO_health_mediation](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md).

- groupings:

  Character vector of groupings to run, any of `"global"`, `"region"`,
  `"income"`, and `"region_income"`. Default
  `c("global", "region", "income")` (the region-by-income cells are the
  slowest and are opt-in).

- error_level:

  Target familywise error rate for the mediator-screening step. Default
  0.05, as in Yu and Kelley (in press).

- lambda_grid:

  Tuning-parameter grid for the subgroup penalized fits (regional,
  income, and region-by-income models). Default
  `seq(0.1, 10, length.out = 100)`, as in Yu and Kelley (in press).

- lambda_grid_global:

  Tuning-parameter grid for the global (all-country) model only. Default
  `seq(0.1, 5, length.out = 100)`, the narrower range Yu and Kelley (in
  press) use for the global fit.

- data:

  The source data. Defaults to
  [WHO_health_mediation](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md).

- cores:

  Number of CPU cores; values above 1 fit the subgroups in parallel via
  the base parallel package (Unix only). Default 1.

- verbose:

  Logical; if `TRUE`, print each group as it is fit. Default `FALSE`.

- full_table:

  Logical; if `TRUE`, return the article's *extended* table layout, with
  the power-enhanced p-value and the active mediator codes reported
  separately under each multiplicity method (Bonferroni for FWER, BH and
  BY for FDR) instead of a single primary method. Default `FALSE`.

## Value

A `data.frame` with one row per fitted subgroup. With
`full_table = FALSE` (the default) the columns are `grouping` (global /
region / income / region_income), `group` (the specific group, for
example `AFR` or `Low`), `n_countries` (number of countries),
`pval_hdmm` and `pval_pe` (the benchmark and power-enhanced p-values for
the total indirect effect), `n_active` (number of active mediators
identified), and `active_mediators` (their codes, or `"none"`). With
`full_table = TRUE` the single `pval_pe`/`n_active`/`active_mediators`
columns are replaced by `pval_pe_bonferroni`, `pval_pe_bh`, `pval_pe_by`
and the matching `active_bonferroni`, `active_bh`, `active_by` code
columns. The result is a `poem_tbl` with the attributes `"outcome"`,
`"full_table"`, and `"notes"` (see Details); call
[summary()](https://yelleknek.github.io/POEM/reference/summary.poem_tbl.md)
on it for a cross-group digest of where the power-enhanced test detects
mediation the benchmark misses.

## Details

The output is the package's **benchmark fit**: applied to the benchmark
[WHO_health_mediation](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md)
data with the default settings, it returns the active mediators reported
in the article, and the package's test suite checks those values, so the
fit serves as a fixed reference point for the method.

Each subgroup is fit with the article's preprocessing: the exposure and
the 57 mediators are standardized, the outcome is centered, and the
confounder indicators (region, income, year) are left unscaled. Groups
with too few complete observations (some region-by-income cells are
empty or nearly so) are skipped and reported as `NA` rows, named in a
message; a group whose fit fails is likewise reported as `NA`, with a
warning that quotes the error. Both are recorded in the `"notes"`
attribute (a `data.frame` with `group`, `stage`, and `message`, or
`NULL` when every group was fit). Because every cell fits a penalized
model over a grid of tuning parameters, a full run with all groupings
fits dozens of models and takes about ten seconds on a current laptop;
start with the default groupings, or a single outcome, before scaling
up.

Exact p-values depend on modeling choices that the article fix (the
tuning-parameter grid, the unscaled confounders); the active mediators
the PE test identifies are the stable, substantive output and reproduce
the article's findings (for example, general government expenditure as a
percent of GDP, `gge_gdp`, for the global infant-mortality model).

## References

Yu, X., and Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

## See also

[`WHO_mediation_design()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_design.md)
for a single design,
[WHO_health_mediation](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md)
for the data,
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
for the test.

Other WHO mediation data:
[`WHO_health_mediation`](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md),
[`WHO_indicator_codebook`](https://yelleknek.github.io/POEM/reference/WHO_indicator_codebook.md),
[`WHO_mediation_design()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_design.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# The global infant-mortality model, the article's headline row (under a
# second). Add "region" and "income" to groupings for the full table.
WHO_mediation_analysis("imr", groupings = "global")
#>  grouping group n_countries pval_hdmm pval_pe  n_active active_mediators
#>  global   ALL   91          0.0409    < 0.0001 1        gge_gdp         
#> 
#> Outcome: imr
```
