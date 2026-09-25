# Broom verbs for POEM results

`tidy()` and `glance()` methods so a POEM result composes with the broom
ecosystem. For a
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
```

## Arguments

- x:

  A `poem_tbl` result.

- ...:

  Unused, for generic compatibility.

## Value

A `data.frame`.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md).

## Author

Xiufan Yu and Ken Kelley
