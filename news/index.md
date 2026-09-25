# Changelog

## POEM 1.0.0

### New features

- [`simulate_guo_mediation()`](https://yelleknek.github.io/POEM/reference/simulate_guo_mediation.md)
  and the `guo_calibration` data object: the manuscript’s
  real-data-motivated heterogeneous setting, with eleven active
  mediators of mixed sign calibrated to the DNA-methylation case study
  of Guo et al. (2022). The total indirect effect is about `-1.597 * c1`
  (the article’s `-1.5977` is from the full-precision coefficients); the
  calibration constants are shipped and inspectable. This is a genuinely
  heterogeneous data-generating process, distinct from the homogeneous
  and contrasting presets of
  [`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md).
- [`pe_identification_study()`](https://yelleknek.github.io/POEM/reference/pe_identification_study.md):
  a Monte Carlo study of how well the power-enhanced screen recovers the
  *individual* active mediators, reporting empirical family-wise error
  rate, false discovery rate, precision, and recall by signal strength
  and multiplicity method (reproduces the manuscript supplement’s
  identification study).
- [`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md)
  now reports the power-enhanced statistic, its enhancement component
  `J_m`, and the global p-value under each multiplicity method
  (Bonferroni / BH / BY) side by side, not just the active sets.
- `WHO_mediation_analysis(full_table = TRUE)` returns the manuscript’s
  extended table layout, with the power-enhanced p-value and active
  mediators reported separately under each multiplicity method.
- A [`summary()`](https://rdrr.io/r/base/summary.html) method for
  `poem_tbl` results: for a grouped comparison table such as the output
  of
  [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md),
  it digests the table across groups and headlines where the
  power-enhanced test detects mediation that the benchmark test misses,
  with the most frequently flagged indicators. Any other `poem_tbl`
  falls back to the ordinary data-frame summary.
- A `drop_constant` argument to
  [`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md)
  and the outcome-specific workers: zero-variance mediator columns are
  dropped with a warning instead of raising an error, and the reported
  active set is renumbered back to the original mediator columns.
  [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md)
  applies this per subgroup, so small region-by-income cells now fit
  instead of being skipped.

### Performance

- The penalized linear fit builds its `[M, X, S]` design matrix once per
  model rather than
  re-[`cbind()`](https://rdrr.io/r/base/cbind.html)-ing it on each of
  the roughly eighty glmnet calls across the lambda grid. The result is
  numerically identical (the published benchmark values are unchanged)
  and about 8% faster per fit.
- The Monte Carlo and empirical functions
  ([`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md),
  [`pe_identification_study()`](https://yelleknek.github.io/POEM/reference/pe_identification_study.md),
  [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md))
  parallelize across replications or subgroups via `cores`, scaling
  near-linearly (about 3.3x on 4 cores). This is the practical lever for
  manuscript-scale studies; the per-fit penalized solver is
  glmnet-bound.
