# Power Enhancement for High-Dimensional Mediation Analysis (POEM)

Powerful global tests for high-dimensional mediation analysis, for
continuous, binary, and count outcomes. The package answers the question
“is there any active mediator among a large set of candidate mediators?”
and, when the answer is yes, reports which individual mediators are
active. The name comes from POwer-Enhanced Mediation.

## The problem POEM solves

A mediation analysis asks how an exposure \\X\\ affects an outcome \\Y\\
through intermediate variables (mediators) \\M\\. With many candidate
mediators, the standard approach tests the *total* indirect effect
\\\beta = \Gamma_x \alpha_m\\, the sum of every mediator's individual
indirect effect. That test has a blind spot: when some mediators carry
positive indirect effects and others carry negative ones, the individual
effects can cancel and the total indirect effect is zero even though
active mediators plainly exist. A test built on the total indirect
effect is then powerless, by construction, against exactly the
heterogeneous mediation that is common in practice.

The power enhancement (PE) tests in this package add a component \\J_m\\
that accumulates the marginal signal from each individual mediator.
Because \\J_m\\ is a sum of magnitudes, opposite-signed indirect effects
reinforce rather than cancel, so the PE test detects active mediators
whether their effects are homogeneous, heterogeneous, or exactly
contrasting, while keeping the Type I error rate of the original test
under the global null of no mediation.

## Main functions

- [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md):

  The front end. Give it the exposure, outcome, mediators, and optional
  confounders, name the outcome type (continuous, binary, or count), and
  it returns a tidy table with the benchmark Wald test, the
  power-enhanced test, the estimated total indirect effect, and the set
  of active mediators it identified.

- [`pe_mediation_linear()`](https://yelleknek.github.io/POEM/reference/pe_mediation_linear.md),
  [`pe_mediation_logistic()`](https://yelleknek.github.io/POEM/reference/pe_mediation_logistic.md),
  [`pe_mediation_poisson()`](https://yelleknek.github.io/POEM/reference/pe_mediation_poisson.md):

  The outcome-specific workers that
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
  dispatches to. Call them directly when the outcome type is known.

- [`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md):

  Generates mediation data under the homogeneous, heterogeneous, and
  contrasting patterns studied in the article, for any of the three
  outcome types.

- [`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md):

  Runs the article's Monte Carlo size and power study: it sweeps a
  signal-strength grid and returns the empirical rejection rate of the
  benchmark and power-enhanced tests at each point.

- [`pe_mediate()`](https://yelleknek.github.io/POEM/reference/pe_mediate.md):

  A data-frame / formula front end to
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
  for users who prefer naming columns to matrices.

- [`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md)
  and
  [`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md):

  The per-mediator screen detail behind a fit, and the active set under
  each multiplicity method.

- [`pe_simulation_study()`](https://yelleknek.github.io/POEM/reference/pe_simulation_study.md)
  and
  [`ss_power_pe_mediation()`](https://yelleknek.github.io/POEM/reference/ss_power_pe_mediation.md):

  Run the size and power study across patterns, and plan the sample size
  for a target power.

- [`plot.poem_tbl()`](https://yelleknek.github.io/POEM/reference/plot.poem_tbl.md),
  [`tidy.poem_tbl()`](https://yelleknek.github.io/POEM/reference/poem_broom.md),
  [`glance.poem_tbl()`](https://yelleknek.github.io/POEM/reference/poem_broom.md):

  A base-graphics plot and broom verbs for POEM results.

- [`mediation_ar1_cov()`](https://yelleknek.github.io/POEM/reference/mediation_ar1_cov.md):

  The autoregressive mediator covariance used throughout the
  simulations, exported as a reusable utility.

- [WHO_health_mediation](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md),
  [`WHO_mediation_design()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_design.md),
  and
  [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md):

  The article's empirical data (economic growth, 57 health-expenditure
  mediators, five health outcomes for 91 WHO members), a helper that
  assembles an analysis design from it, and a function that reproduces
  the article's full set of subgroup analyses.

## Reading the output

Every estimation and testing function returns a tidy `data.frame` with a
`term` column and a numeric `value` column, carrying the `poem_tbl`
class so it prints with whole numbers shown without a decimal part,
other quantities to a few significant figures, and p-values to four
decimal places. The stored numbers keep full precision; only the display
rounds. See
[poem_tbl](https://yelleknek.github.io/POEM/reference/poem_tbl.md) for
the details.

## References

Yu, X., and Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

Fan, J., Liao, Y., and Yao, J. (2015). Power enhancement in
high-dimensional cross-sectional tests. *Econometrica, 83*(4),
1497–1541. [doi:10.3982/ECTA12749](https://doi.org/10.3982/ECTA12749)

Guo, X., Li, R., Liu, J., and Zeng, M. (2022). High-dimensional
mediation analysis for selecting DNA methylation loci mediating
childhood trauma and cortisol stress reactivity. *Journal of the
American Statistical Association, 117*(539), 1110–1121.
[doi:10.1080/01621459.2022.2053136](https://doi.org/10.1080/01621459.2022.2053136)

## See also

Useful links:

- <https://github.com/yelleKneK/POEM>

- <https://yelleknek.github.io/POEM/>

- Report bugs at <https://github.com/yelleKneK/POEM/issues>

## Author

Xiufan Yu and Ken Kelley
