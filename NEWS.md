# POEM 1.0.0

## Tuning, empty fits, and the Monte Carlo studies

* The default tuning grid for a continuous outcome is now the grid the
  article's simulations used, rescaled to the `n` and `p` at hand by the
  rate `sqrt(log(p) / n)` the theory requires of the SCAD tuning
  parameter (`pe_lambda_grid()`, new). The review-era code searched a
  fixed grid of 20 values from 0.05 to 1 for every design; at the
  article's linear design that grid lets the identified set's familywise
  error rate rise from the reported 0.000 to about 0.17, and the article's
  linear tables did not reproduce under it. They do now, within Monte
  Carlo error. The price is global-test power at small designs, measured
  and stated on `?pe_lambda_grid`; pass the review-era grid explicitly to
  trade the guarantee back. Binary and count outcomes keep the review-era
  grid as their default until the article's per-panel grids for those
  families are on record.
* Every fit records the grid it searched and the value HBIC chose in a
  `"tuning"` attribute and reports it in the print footer, with a note
  when the choice sat at the grid's lower end (where the grid's floor,
  not HBIC, decided the fit).
* A penalized fit that selects no mediator at any grid value now warns
  (condition class `poem_empty_fit`), sets an `"empty_fit"` attribute,
  says so in the print footer, and keeps the full row schema with `NA`
  confidence limits, instead of silently returning statistics of 0 and
  p-values of 1 with no interval rows. The values follow the article's
  Monte Carlo convention (an empty selection is a non-rejection).
* `pe_power_curve()`, `pe_simulation_study()`, `pe_identification_study()`,
  and `ss_power_pe_mediation()` accept `lambda_grid` and
  `lambda_grid_reduced`, so the article's studies can be run at the grids
  behind its tables; they report the empty fits per grid point in a new
  `n_empty` column, validate their design arguments instead of returning
  `NaN` rates, warn when replications fail (quoting the first error), and
  refuse an `outcome_args` entry (such as `seed`) that would override the
  design or make every replication identical.
* `pe_identification_study()` scores the four metrics as the article's
  supplement does (a mediator with a nonzero outcome coefficient counts as
  active, and an empty selection scores precision 0), with a `truth`
  argument for the both-paths definition. The supplement's tables now
  reproduce within Monte Carlo error under the default grid.
* The confidence-interval rows `total_indirect_lower` and
  `total_indirect_upper` (a Wald interval) are documented, and the
  interval rows are present on every fit.
* Seeds go through `withr::local_seed()`: a supplied `seed` reproduces a
  run bit for bit and the caller's generator state and kind are restored
  on exit, with `cores > 1` included; the package never touches
  `.Random.seed` or the global environment itself.
* Input validation names the argument: missing, `NaN`, or infinite values
  in `X`, `Y`, `M`, or `Z`, a constant exposure or confounder column, and
  a single-column `M` stop with a clear message instead of surfacing as an
  error from glmnet or the coordinate descent.
* The broom generics `tidy()` and `glance()` are re-exported, so they work
  after `library(POEM)` alone; `glance()` reports one total-indirect-effect
  column per exposure. `pe_selection()` returns a `poem_tbl`, and the
  p-value columns of every wide table (`pval_hdmm`, `pval_pe`, the
  per-method columns, `screen_p`) print to four fixed decimals with the
  `< 0.0001` floor like the p-value rows of a fit.
* The test suite pins the three shipped example fits to the values of the
  reference implementation at full precision, pins the article's
  real-data active sets exactly, and verifies the seed contract of every
  seeded function at runtime.
* Every reference carries its DOI, and the article's sections are cited
  by the numbering of the accepted article.

## Names, arguments, and vocabulary

* The level at which a replication or a group counts as a rejection is
  `alpha_level` in `pe_power_curve()`, `pe_simulation_study()`,
  `ss_power_pe_mediation()`, and `summary.poem_tbl()` (it was `alpha`),
  the portfolio's name for a significance level; `error_level` remains the
  screening level of the mediator-identification step.
* Column names travel to the output. With a named `M` (a named matrix or a
  data.frame) the print footer names the active mediators, the mediator
  table (`pe_mediators()`, `tidy()`) gains a `name` column, and
  `pe_selection()` lists names; the `"mediator_names"` attribute carries
  them. With several named exposures the per-exposure rows are suffixed
  with the exposure's column name instead of `_1`, `_2`, ... . An unnamed
  analysis is unchanged, and no number or position changes.
* `simulate_mediation_data()` validates `n`, `p`, `q`, `d`, `c1`, `c2`,
  `rho`, and `sigma_y` by name; `WHO_mediation_analysis()` validates its
  levels, grids, `cores`, flags, and `data` up front (a wrong `data` used
  to return a table of `NA` rows without a word), names the groups it
  skipped in a message and any failed fit in a warning, and records both
  in a `"notes"` attribute; `summary.poem_tbl()` validates `alpha_level`;
  `print.poem_tbl()` and `format.poem_tbl()` validate `digits` and
  `digits_p`.
* The documentation says "the article" and "Yu and Kelley (in press)"
  throughout (the paper is in press, no longer a manuscript), the
  reproduction vignette is `reproducing-the-article`, the name expands one
  way (Power Enhancement for High-Dimensional Mediation Analysis; from
  POwer-Enhanced Mediation), and the error rate is "familywise".
* `?pe_mediation` and the guided tour state when POEM is the right tool
  and when it is not (small `n`, strongly correlated mediators, rare
  binary outcomes, overdispersed counts, missing data), and README points
  to DMAR for low-dimensional mediation models. Every vignette closes with
  its session information.

## New features

* `simulate_guo_mediation()` and the `guo_calibration` data object: the
  article's real-data-motivated heterogeneous setting, with eleven
  active mediators of mixed sign calibrated to the DNA-methylation case
  study of Guo et al. (2022). The total indirect effect is about
  `-1.597 * c1` (the article's `-1.5977` is from the full-precision
  coefficients); the calibration constants are shipped and inspectable.
  This is a genuinely heterogeneous data-generating process, distinct from
  the homogeneous and contrasting presets of `simulate_mediation_data()`.
* `pe_identification_study()`: a Monte Carlo study of how well the
  power-enhanced screen recovers the *individual* active mediators, reporting
  empirical familywise error rate, false discovery rate, precision, and
  recall by signal strength and multiplicity method (reproduces the
  article supplement's identification study).
* `pe_selection()` now reports the power-enhanced statistic, its enhancement
  component `J_m`, and the global p-value under each multiplicity method
  (Bonferroni / BH / BY) side by side, not just the active sets.
* `WHO_mediation_analysis(full_table = TRUE)` returns the article's
  extended table layout, with the power-enhanced p-value and active
  mediators reported separately under each multiplicity method.
* A `summary()` method for `poem_tbl` results: for a grouped comparison
  table such as the output of `WHO_mediation_analysis()`, it digests the
  table across groups and headlines where the power-enhanced test detects
  mediation that the benchmark test misses, with the most frequently
  flagged indicators. Any other `poem_tbl` falls back to the ordinary
  data-frame summary.
* A `drop_constant` argument to `pe_mediation()` and the outcome-specific
  workers: zero-variance mediator columns are dropped with a warning
  instead of raising an error, and the reported active set is renumbered
  back to the original mediator columns. `WHO_mediation_analysis()` applies
  this per subgroup, so small region-by-income cells now fit instead of
  being skipped.

## Performance

* The penalized linear fit builds its `[M, X, S]` design matrix once per
  model rather than re-`cbind()`-ing it on each of the roughly eighty
  glmnet calls across the lambda grid. The result is numerically identical
  (the published benchmark values are unchanged) and about 8% faster per
  fit.
* The Monte Carlo and empirical functions (`pe_power_curve()`,
  `pe_identification_study()`, `WHO_mediation_analysis()`) parallelize
  across replications or subgroups via `cores`, scaling near-linearly
  (about 3.3x on 4 cores). This is the practical lever for
  article-scale studies; the per-fit penalized solver is glmnet-bound.
