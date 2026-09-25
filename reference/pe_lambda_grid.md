# Default tuning-parameter grid for the penalized mediator fit

Builds the grid of candidate SCAD tuning parameters that
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
searches by the high-dimensional BIC when no `lambda_grid` is supplied.
The grid is the one the article's simulation scripts used at their
design (\\n = 300\\, \\p = 500\\), rescaled to the sample size and
number of mediators at hand by the rate the theory requires of the
tuning parameter, \\\sqrt{\log p / n}\\.

## Usage

``` r
pe_lambda_grid(
  n,
  p,
  outcome = c("continuous", "binary", "count"),
  model = c("full", "reduced"),
  length.out = NULL
)
```

## Arguments

- n, p:

  Number of observations and candidate mediators.

- outcome:

  Outcome type: `"continuous"`, `"binary"`, or `"count"`.

- model:

  Which fit the grid is for: `"full"` (the mediator model used for
  selection and the power-enhanced screen) or `"reduced"` (the reduced
  model the continuous-outcome Wald test refits; for binary and count
  outcomes the reduced grid is not used and `"reduced"` returns the same
  grid as `"full"`).

- length.out:

  Number of grid values. Default `NULL` uses the article's counts: 20
  for a continuous outcome, 15 for binary and count outcomes.

## Value

A numeric vector of increasing candidate tuning parameters.

## Details

The SCAD penalty in Yu and Kelley (in press) has a single tuning
parameter \\\lambda\\; the fit is repeated at each value of a candidate
grid and the value minimizing the high-dimensional BIC (HBIC) is kept.
Two facts make the grid, and especially its lower end, part of the
estimator rather than a detail. First, the article's theory requires
\\\lambda\\ to be of larger order than \\\sqrt{\log p / n}\\ (and than
\\\sqrt{s_0 / n}\\, with \\s_0\\ the number of active mediators); a grid
whose values fall below that regime lets HBIC admit spurious mediators
to the penalized selection, and the screening step then reports some of
them as active. Second, the HBIC minimum drifts toward small \\\lambda\\
and often sits at the smallest value offered, so the grid's lower end
frequently decides the fit.

The article's scripts set the grid by hand for each outcome family at
\\n = 300\\, \\p = 500\\: for the continuous outcome, 20 values from 0.2
to 0.39 (full model) and from 0.27 to 0.46 (reduced model); for the
binary outcome, 15 values from 0.04 to 0.2; for the count outcome, 15
values from 0.7 to 5. At that design \\\sqrt{\log p / n} = 0.1439\\, so
the continuous grid's lower end is about 1.4 times the rate. This
function returns exactly those grids at the article's design and scales
their endpoints by \\\sqrt{\log p / n}\\ elsewhere, which keeps the
lower end in the regime the theory requires as \\n\\ and \\p\\ change.
The multipliers differ across families because the penalized likelihood
of a logistic or Poisson outcome is on a different scale from penalized
least squares.

Which grid
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
searches by default depends on the outcome. For a continuous outcome the
default is this grid. The review-era implementation searched a fixed
grid of 20 values from 0.05 to 1 for every family and design; at the
article's linear design that grid lets the identified set's familywise
error rate rise from the reported 0.000 to about 0.17 under contrasting
mediation (0.35 to 0.39 with a finer grid over the same low range), and
the article's linear tables do not reproduce under it. The guarantee has
a price at small designs, where the rescaled floor removes more signal:
at \\n = 100\\, \\p = 40\\ under contrasting mediation with \\c_1 = 1\\
the power-enhanced test rejects in 17 percent of replications with this
grid against 77 percent with the review-era grid, while the identified
set's familywise error rate is 0.000 against 0.07; at \\n = 200\\, \\p =
100\\ the figures are 58 against 97 percent and 0.000 against 0.13 (200
replications each). Pass `lambda_grid = seq(0.05, 1, length.out = 20)`
to trade the guarantee for that power, or to reproduce results computed
with the review-era grid; the print footer reports which value HBIC
chose.

For binary and count outcomes the default of
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
remains the review-era grid. The article's scripts for those families
carry one grid each (0.04 to 0.2 and 0.7 to 5), which reproduce the
article's homogeneous panels, while its contrasting panels reproduce
under the review-era grid and not under the family grids (the Poisson
grid leaves about a third of contrasting-count fits with no selected
mediator); the grids behind each panel were set by hand and not
recorded. This function still returns the family scripts' grids for
those outcomes, for use with the homogeneous designs, and the default
will follow the article's grids once they are on record.

## References

Yu, X., and Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

Wang, L., Kim, Y., and Li, R. (2013). Calibrating nonconvex penalized
regression in ultra-high dimension. *The Annals of Statistics, 41*(5),
2505–2536. [doi:10.1214/13-AOS1159](https://doi.org/10.1214/13-AOS1159)

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
whose `lambda_grid` and `lambda_grid_reduced` arguments default to this
grid.

Other mediation tests:
[`pe_mediate()`](https://yelleknek.github.io/POEM/reference/pe_mediate.md),
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
# The article's linear-model grid at its own design
pe_lambda_grid(n = 300, p = 500, outcome = "continuous")
#>  [1] 0.20 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.30 0.31 0.32 0.33 0.34
#> [16] 0.35 0.36 0.37 0.38 0.39
# The same grid rescaled to a smaller study
pe_lambda_grid(n = 100, p = 200, outcome = "continuous")
#>  [1] 0.3198547 0.3358474 0.3518402 0.3678329 0.3838256 0.3998184 0.4158111
#>  [8] 0.4318038 0.4477966 0.4637893 0.4797820 0.4957748 0.5117675 0.5277603
#> [15] 0.5437530 0.5597457 0.5757385 0.5917312 0.6077239 0.6237167
pe_lambda_grid(n = 300, p = 500, outcome = "count")
#>  [1] 0.700000 1.007143 1.314286 1.621429 1.928571 2.235714 2.542857 2.850000
#>  [9] 3.157143 3.464286 3.771429 4.078571 4.385714 4.692857 5.000000
```
