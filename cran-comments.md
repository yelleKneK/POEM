# cran-comments

## Submission

This is a new submission: POEM 1.0.0.

POEM implements the methods of the accompanying peer-reviewed article,
Yu, X., & Kelley, K. (in press), "Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis," *Journal of the American Statistical
Association*, whose abstract names the software POEM.

## Package name

`R CMD check --as-cran` notes that the name POEM is similar, case
insensitively, to the Bioconductor package `poem` ("POpulation-based
Evaluation Metrics"). We respectfully request the name POEM:

* POEM is the name the JASA article gives this software, so it is already
  established in the literature. The name comes from POwer-Enhanced
  Mediation.
* POEM is not on CRAN, and the Bioconductor `poem` is in an unrelated
  domain (evaluation metrics for population-level analyses, not
  mediation analysis or hypothesis testing); the two names also differ in
  case.

If CRAN nonetheless prefers a different name, we are happy to discuss
alternatives.

## Test environments

* local macOS (Apple Silicon), R 4.6.1 (2026-09-25): full
  `R CMD check --as-cran` with both manuals built, from the release
  tarball, built from a clean archive of the repository by the
  maintainer's release gate.
* win-builder R-devel and R-release: to be run on this tarball before
  upload, with the date and status recorded here.
* GitHub Actions R-CMD-check matrix (Ubuntu release and devel, macOS,
  Windows, Ubuntu oldrel), `.github/workflows/R-CMD-check.yaml`, on the
  pushed commit.

## R CMD check results

0 errors, 0 warnings. The expected NOTEs:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Ken Kelley <kkelley@nd.edu>'

New submission

Conflicting package names (submitted: POEM, existing: poem [Bioconductor])
```

The "New submission" note is expected for a first submission; the name
conflict is addressed above.

On the local machine the check also reports "Skipping checking HTML
validation: 'tidy' doesn't look like recent enough HTML Tidy"; that is
environmental (the CRAN check farm has a current HTML Tidy).

## Check time

The package contains no `\donttest{}` and no `\dontrun{}`, so the
examples run in a single pass: 9 s in the check across the 24 help pages
with examples, the slowest page (`summary.poem_tbl`) about 1 s. The
tests run in 19 s on CRAN's path (the Monte Carlo blocks and the full
real-data reproduction carry `skip_on_cran()` and run in the
maintainer's release gate, where the whole suite of 440 expectations
passes) and the three vignettes build in 105 s. The whole check took
167 s locally (2026-09-25), the feasibility step excluded.

## Possibly misspelled words in DESCRIPTION

* **POwer** is deliberate: it shows where the package name's letters come
  from (POwer-Enhanced Mediation).

## Downstream dependencies

This is a new package; no reverse dependencies on CRAN.
