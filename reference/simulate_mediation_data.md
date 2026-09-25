# Simulate high-dimensional mediation data

Generates a data set from the mediation models studied in the article,
for a continuous, binary, or count outcome, under the homogeneous or
contrasting (heterogeneous) mediation patterns. The mediators carry an
AR(1) correlation structure, a handful are truly active, and the rest
are null, so the data exercise exactly the situation the power-enhanced
tests are designed for. This is the generator behind the package's
example data sets and behind
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md).

## Usage

``` r
simulate_mediation_data(
  n,
  p,
  outcome = c("continuous", "binary", "count"),
  pattern = c("homogeneous", "contrasting", "heterogeneous"),
  c1 = 1,
  c2 = 0.5,
  q = 1L,
  d = 0L,
  alpha_m = NULL,
  tau = NULL,
  rho = 0.5,
  sigma_y = 0.5,
  alpha_z = NULL,
  Gamma_z = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Number of observations.

- p:

  Number of candidate mediators.

- outcome:

  Outcome type: `"continuous"`, `"binary"`, or `"count"`.

- pattern:

  Mediation pattern: `"homogeneous"` (active effects share a sign) or
  `"contrasting"` (active effects cancel; total indirect effect zero).
  `"heterogeneous"` is accepted as a synonym for `"contrasting"`.
  Ignored when `alpha_m` is supplied directly.

- c1:

  Scale on the exposure-on-mediator coefficients \\\Gamma_x = c_1
  \tau\\. `c1 = 0` makes every mediator inactive (the global null).
  Default 1.

- c2:

  The direct effect \\\alpha_x\\ of each exposure on the outcome.
  Default 0.5.

- q:

  Number of exposures. Default 1. For `q > 1` you must supply `tau` as a
  \\q \times p\\ matrix.

- d:

  Number of confounders. Default 0 (none). When positive, confounders
  are drawn standard normal and enter both models through `alpha_z` and
  `Gamma_z`.

- alpha_m:

  Optional length-\\p\\ vector of mediator-on-outcome coefficients,
  overriding the `pattern` preset.

- tau:

  Optional base pattern for \\\Gamma_x\\. Default builds the article's
  \\\tau = (0.1, 0.2, 0.3, 0.4, 0.5, \text{noise})\\ for `q = 1`;
  required (as a \\q \times p\\ matrix) when `q > 1`.

- rho:

  AR(1) correlation parameter for the mediator noise. Default 0.5.

- sigma_y:

  Outcome noise standard deviation for a continuous outcome. Default
  0.5.

- alpha_z, Gamma_z:

  Optional confounder coefficients (length \\d\\ vector and \\d \times
  p\\ matrix). Default zero when `d > 0`.

- seed:

  Optional integer seed. When supplied it is set locally and the
  caller's random number generator state is restored on exit; `NULL`
  (default) leaves the random number generator alone.

## Value

A list with the simulated data and the truth used to generate it: `X`
(\\n \times q\\), `M` (\\n \times p\\), `Y` (length \\n\\), `Z` (\\n
\times d\\ or `NULL`), the coefficient values `alpha_m`, `Gamma_x`,
`alpha_x`, the total indirect effect `beta = Gamma_x alpha_m`, and
`active_mediators` (the indices that are truly active, i.e. nonzero in
both paths).

## Details

One exposure column (or \\q\\ of them) is drawn standard normal. The
mediators follow \\M = X \Gamma_x + Z \Gamma_z + \varepsilon_m\\ with
\\\varepsilon_m \sim N(0, \Sigma)\\, \\\Sigma\\ the AR(1) covariance
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEM/reference/mediation_ar1_cov.md)
with parameter `rho`. The exposure-on-mediator coefficient is \\\Gamma_x
= c_1 \tau\\, where the base pattern \\\tau\\ has small increasing
loadings on the first five mediators and random noise loadings on the
rest, so `c1` scales the overall \\X \to M\\ signal (and `c1 = 0` gives
the global null of no mediation). The outcome is generated from its
model with mediator coefficients \\\alpha_m\\ (the preset for the chosen
`pattern`, or a user-supplied vector) and direct effect \\\alpha_x =
c_2\\:

- continuous: \\Y = M\alpha_m + X\alpha_x + Z\alpha_z + \varepsilon_y\\,
  \\\varepsilon_y \sim N(0, \sigma_y^2)\\;

- binary: \\Y \sim \mathrm{Bernoulli}(\mathrm{logit}^{-1}( M\alpha_m +
  X\alpha_x + Z\alpha_z))\\;

- count: \\Y \sim \mathrm{Poisson}(\exp(\eta))\\, with the log-mean
  \\\eta\\ clamped to \\\[-5, 5\]\\ to guard against overflow, as in the
  article.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
to analyze the data,
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
to sweep `c1`,
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEM/reference/mediation_ar1_cov.md)
for the mediator covariance.

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEM/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEM/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEM/reference/pe_identification_study.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md),
[`pe_simulation_study()`](https://yelleknek.github.io/POEM/reference/pe_simulation_study.md),
[`plot.poem_tbl()`](https://yelleknek.github.io/POEM/reference/plot.poem_tbl.md),
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEM/reference/simulate_guo_mediation.md),
[`ss_power_pe_mediation()`](https://yelleknek.github.io/POEM/reference/ss_power_pe_mediation.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 100, p = 50, outcome = "continuous",
                             pattern = "contrasting")
dim(d$M)
#> [1] 100  50
d$active_mediators       # truly active mediators
#> [1] 1 2 3 4
d$beta                   # total indirect effect (zero under "contrasting")
#> [1] 0
```
