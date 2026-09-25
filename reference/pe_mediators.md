# Per-mediator detail behind a power-enhanced mediation test

Returns the per-mediator table that explains a
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
result: one row for each mediator the penalized fit selected as a
candidate, with its two path statistics, the screening p-value, and
whether it was identified as active. This is the "why" behind the active
set, so a user can see which mediators were close to the screening
threshold and which were screened out.

## Usage

``` r
pe_mediators(fit)
```

## Arguments

- fit:

  A result from
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
  or one of its workers.

## Value

A tidy `data.frame` (class `poem_tbl`) with columns:

- mediator:

  Column position of the mediator in the original `M`.

- name:

  The mediator's column name in `M`; present only when `M` has column
  names.

- t_outcome:

  The standardized mediator-on-outcome statistic (the \\M \to Y\\ path,
  \\\hat\alpha\_{m,j} / \hat\sigma\_{m,j}\\).

- t_exposure:

  The strongest standardized exposure-on-mediator statistic across
  exposures (the \\X \to M\\ path).

- screen_p:

  The screening p-value, the larger of the two path p-values (the
  smaller across exposures when there are several). A mediator is active
  when this clears the multiplicity threshold.

- selected:

  Logical; whether the mediator was identified as active.

The table is empty when the fit selected no candidate mediators.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md).

Other mediation tests:
[`pe_lambda_grid()`](https://yelleknek.github.io/POEM/reference/pe_lambda_grid.md),
[`pe_mediate()`](https://yelleknek.github.io/POEM/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md),
[`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md),
[`summary.poem_tbl()`](https://yelleknek.github.io/POEM/reference/summary.poem_tbl.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                             pattern = "contrasting", c1 = 1)
fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
pe_mediators(fit)
#>  mediator t_outcome t_exposure screen_p selected
#>  1        18.22     0.3641     0.7158   FALSE   
#>  2        -7.179    3.131      0.0017    TRUE   
```
