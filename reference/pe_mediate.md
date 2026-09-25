# Power-enhanced mediation test from a data frame (formula interface)

A convenience front end to
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
for users who have a data frame rather than ready-made matrices. Give it
the outcome and exposure as a formula, name the mediator columns, and
(optionally) name confounder columns; it assembles the exposure,
outcome, mediator, and confounder matrices (turning factor confounders
into indicator variables) and runs the test. This is the gentler entry
point; the matrix interface
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
is the workhorse.

## Usage

``` r
pe_mediate(
  formula,
  data,
  mediators,
  confounders = NULL,
  outcome = c("continuous", "binary", "count"),
  ...
)
```

## Arguments

- formula:

  A two-sided formula `outcome ~ exposure` naming the outcome and one or
  more exposure columns in `data`.

- data:

  A `data.frame` containing the outcome, exposure, mediator, and any
  confounder columns.

- mediators:

  Character vector of mediator column names in `data` (the candidate
  mediators `M`).

- confounders:

  Optional character vector of confounder column names. Factor and
  character columns are expanded into indicator variables; numeric
  columns enter as is. Default `NULL`.

- outcome:

  Outcome type: `"continuous"`, `"binary"`, or `"count"`.

- ...:

  Further arguments passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
  (for example `method`, `error_level`, `conf_level`,
  `report_all_methods`, `lambda_grid`).

## Value

The tidy `poem_tbl` returned by
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md).

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
for the matrix interface,
[`WHO_mediation_design()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_design.md)
for the worked health-expenditure design.

Other mediation tests:
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
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
# Build a small data frame and test through the formula interface.
set.seed(113)
d <- simulate_mediation_data(n = 150, p = 40, outcome = "continuous",
                             pattern = "contrasting", c1 = 1)
df <- data.frame(y = d$Y, x = d$X[, 1], d$M)
meds <- grep("^X", names(df), value = TRUE)   # the mediator columns
pe_mediate(y ~ x, data = df, mediators = meds, outcome = "continuous")
#>  term                  value   
#>  stat_hdmm             0.1615  
#>  pval_hdmm             0.6878  
#>  stat_pe               440.3   
#>  j_pe                  440.1   
#>  pval_pe               < 0.0001
#>  total_indirect_effect -0.03065
#>  total_indirect_lower  -0.1801 
#>  total_indirect_upper  0.1188  
#>  n_active_mediators    2       
#>  df                    1       
#>  n_candidate_mediators 40      
#>  n_observations        150     
#> 
#> Outcome model: continuous (linear)
#> Active mediators identified (2): 3, 4
```
