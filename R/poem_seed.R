# Seed discipline for every function that takes a `seed` argument.
#
# A supplied seed is set for the duration of the calling function only. The
# caller's random number generator state, including the generator kind, is
# restored when that function exits, so a call made reproducible by `seed`
# does not perturb the random stream of the script around it. The saving and
# restoring is delegated to withr::local_seed(), the accepted way to do this
# on CRAN (the review that returned the DMAR package on 2026-09-04 named the
# hand-rolled assign-to-globalenv idiom as modifying the global environment);
# the package itself never reads or writes `.Random.seed`. A NULL seed is a
# no-op: the draws come from the user's current generator state, exactly as
# if the function had no `seed` argument.
#
# `parallel = TRUE` runs the call on the L'Ecuyer-CMRG generator, whose
# streams parallel::mclapply() hands out to forked replications so a seeded
# parallel run is reproducible at the same core count. With a seed the
# generator kind is set and restored together with the state by withr; with
# no seed only the kind is switched, and switched back on exit, so the user's
# stream keeps advancing through the call exactly as it did before (an
# unseeded parallel call in a seeded wrapper such as pe_simulation_study()
# must not hand every pattern the same draws).
#
# Call it from the body of the function whose exit should restore the state
# (the default `envir` is that function's frame), never from a helper.
#' @keywords internal
#' @noRd
.poem_local_seed <- function(seed, parallel = FALSE, envir = parent.frame()) {
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))
      stop("`seed` must be a single number or NULL.", call. = FALSE)
    withr::local_seed(as.integer(seed), .local_envir = envir,
                      .rng_kind = if (isTRUE(parallel)) "L'Ecuyer-CMRG" else NULL)
  } else if (isTRUE(parallel)) {
    old_kind <- RNGkind("L'Ecuyer-CMRG")
    withr::defer(do.call(RNGkind, as.list(old_kind)), envir = envir)
  }
  invisible(seed)
}
