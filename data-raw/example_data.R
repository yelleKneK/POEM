# Regenerate the three example data sets shipped with POEM. Each is produced by
# the package's own simulator under a homogeneous mediation pattern with a few
# confounders, at a fixed seed, so the data are reproducible and the "truth"
# (which mediators are active, the total indirect effect) travels with them.
# Run with: source("data-raw/example_data.R") from the package root, then
# document the objects in R/data_examples.R and rebuild.
#
# The datasets are kept deliberately small (n a few hundred, p = 100) so the
# help-page examples and tests that use them run in a few seconds.

devtools::load_all(".")

# Continuous outcome, homogeneous (same-sign) mediation, three confounders.
set.seed(113)
example_continuous <- simulate_mediation_data(
  n = 200, p = 100, outcome = "continuous", pattern = "homogeneous",
  c1 = 0.6, c2 = 0.5, d = 3)

# Binary outcome, homogeneous mediation. A larger sample and signal, since a
# binary outcome carries less information per observation.
set.seed(113)
example_binary <- simulate_mediation_data(
  n = 300, p = 100, outcome = "binary", pattern = "homogeneous",
  c1 = 0.5, c2 = 1, d = 3)

# Count outcome, homogeneous mediation.
set.seed(113)
example_count <- simulate_mediation_data(
  n = 300, p = 100, outcome = "count", pattern = "homogeneous",
  c1 = 0.4, c2 = 0.4, d = 3)

# Verify each truly has active mediators and the expected structure before
# saving, so a regeneration that silently produced a null data set is caught.
stopifnot(
  length(example_continuous$active_mediators) >= 1L,
  length(example_binary$active_mediators)     >= 1L,
  length(example_count$active_mediators)      >= 1L,
  ncol(example_continuous$M) == 100L,
  all(example_binary$Y %in% c(0, 1)),
  all(example_count$Y >= 0 & example_count$Y == round(example_count$Y)))

usethis_available <- requireNamespace("usethis", quietly = TRUE)
if (usethis_available) {
  usethis::use_data(example_continuous, overwrite = TRUE, compress = "xz")
  usethis::use_data(example_binary, overwrite = TRUE, compress = "xz")
  usethis::use_data(example_count, overwrite = TRUE, compress = "xz")
} else {
  save(example_continuous, file = "data/example_continuous.rda", compress = "xz")
  save(example_binary,     file = "data/example_binary.rda",     compress = "xz")
  save(example_count,      file = "data/example_count.rda",      compress = "xz")
}
