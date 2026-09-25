# Tidy printing for POEM result tables

Every estimation and testing function in POEM returns a tidy
`data.frame` with a `term` column and one or more numeric columns. A
single numeric column routinely holds quantities on very different
scales: a whole-number count of mediators next to a chi-square statistic
in the hundreds next to a p-value of \\10^{-12}\\. The base
[`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html) method
formats a whole column with one common format, which forces either a
wall of trailing zeros or a slide into scientific notation. The
`poem_tbl` class supplies `print` and `format` methods that format each
value on its own terms: whole numbers (counts, sample sizes, mediator
indices) print without a decimal part, other values print to a few
significant figures, and p-values print to a fixed number of decimal
places with a “\< 0.0001” floor.

## Usage

``` r
# S3 method for class 'poem_tbl'
format(x, digits = getOption("poem.digits", 4L), digits_p = 4L, ...)

# S3 method for class 'poem_tbl'
print(x, digits = getOption("poem.digits", 4L), digits_p = 4L, ...)
```

## Arguments

- x:

  A `poem_tbl` object (a tidy `data.frame` returned by a POEM function).

- digits:

  Number of significant figures for non-integer values, a single
  positive whole number. Defaults to `getOption("poem.digits", 4L)`.

- digits_p:

  Number of decimal places for p-values, a single positive whole number.
  Defaults to 4. A p-value below `10^(-digits_p)` prints as “\< 0.0001”.

- ...:

  Additional arguments passed to
  [`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html).

## Value

`print.poem_tbl` returns `x` invisibly. `format.poem_tbl` returns a
`data.frame` whose numeric columns have been formatted to character for
display.

## Details

The stored numeric values are never rounded. Only the display changes,
so any arithmetic you do on the returned object (a confidence interval
width, a further calculation) uses full precision. To see more digits,
raise `digits` or read the column directly with `x$value`.

This mirrors the display convention of the DMAR package, whose house
style POEM follows.

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
                             outcome = "continuous")
fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
fit                       # rounded for reading, with the p-value floor
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
fit$value[fit$term == "pval_pe"]   # full precision underneath
#> [1] 9.154634e-40
```
