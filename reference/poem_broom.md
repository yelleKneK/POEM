# Broom verbs for POEM results

`tidy()` and `glance()` methods so a POEM result composes with the broom
ecosystem; the generics are re-exported, so `tidy(fit)` and
`glance(fit)` work after
[`library(POEM)`](https://github.com/yelleKneK/POEM) alone. For a
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
fit, `tidy()` returns the per-mediator table (see
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md))
and `glance()` returns a one-row summary of the two tests and the total
indirect effect. For a
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
table both return the table itself (it is already one row per grid
point).

## Usage

``` r
# S3 method for class 'poem_tbl'
tidy(x, ...)

# S3 method for class 'poem_tbl'
glance(x, ...)

tidy(x, ...)

glance(x, ...)
```

## Arguments

- x:

  A `poem_tbl` result.

- ...:

  Unused, for generic compatibility.

## Value

A `data.frame`. For a fit, `tidy()` has the columns of
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md)
and `glance()` the columns `stat_hdmm`, `pval_hdmm`, `stat_pe`,
`pval_pe`, `total_indirect_effect` (one column per exposure when there
are several, suffixed with the exposure's column name or, for an unnamed
`X`, `_1`, `_2`, ...), `n_active_mediators`, and `n_observations`.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md).

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 120, p = 40, pattern = "contrasting",
                             outcome = "continuous", c1 = 1)
fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
tidy(fit)
#>   mediator t_outcome t_exposure   screen_p selected
#> 1        1  11.19411   2.431541 0.01503477     TRUE
glance(fit)
#>   stat_hdmm  pval_hdmm  stat_pe      pval_pe total_indirect_effect
#> 1  3.362122 0.06671204 175.5098 4.633246e-40             0.1150131
#>   n_active_mediators n_observations
#> 1                  1            120
```
