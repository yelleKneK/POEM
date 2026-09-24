# Internal utilities shared by the POEM testing and simulation functions.
# None are exported. They concentrate the argument coercion, validation, and
# preprocessing that every pe_mediation_*() worker needs, so the workers
# differ only in the parts that are genuinely family-specific (the penalized
# fit and the Wald inference). The prior PEmediation package repeated this
# logic in each worker; centralizing it here removes that duplication and
# gives one place to test the input contract.

# Coerce X / M / Z to a numeric double matrix. A bare numeric vector becomes a
# one-column matrix (the q = 1 single-exposure case, which is the common one),
# a data.frame is matricized, and storage mode is forced to double so the
# downstream linear algebra never silently promotes from integer. NULL passes
# through as NULL so an absent confounder matrix Z stays absent.
#' @keywords internal
#' @noRd
.as_numeric_matrix <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (is.vector(x) && !is.list(x)) x <- matrix(x, ncol = 1L)
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sprintf("`%s` must be a numeric matrix, data.frame, or numeric vector.",
                 name), call. = FALSE)
  x <- as.matrix(x)
  if (!is.numeric(x))
    stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
  storage.mode(x) <- "double"
  x
}

# Coerce the outcome Y to a numeric vector. A one-column matrix is flattened
# (a convenience, since Y is naturally a column), anything list-like is
# rejected, and the result is cast to double.
#' @keywords internal
#' @noRd
.as_numeric_vector <- function(x, name) {
  if (is.matrix(x) && ncol(x) == 1L) x <- as.vector(x)
  if (!is.atomic(x) || is.list(x))
    stop(sprintf("`%s` must be a numeric vector.", name), call. = FALSE)
  as.numeric(x)
}

# The tuning-parameter grid must be a non-empty vector of strictly positive,
# finite values: each entry is a candidate lambda for the penalized fit, and a
# nonpositive or non-finite lambda is not a valid penalty level.
#' @keywords internal
#' @noRd
.validate_lambda_grid <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x)) || any(x <= 0))
    stop(sprintf("`%s` must be a non-empty numeric vector with all elements > 0.",
                 name), call. = FALSE)
  invisible(TRUE)
}

# Resolve the tuning grids a worker will search. With no user grid, a
# continuous outcome searches the article's grids rescaled to this design
# (pe_lambda_grid()); a binary or count outcome searches the review-era grid,
# because the article's per-panel grids for those families are not on record
# (see ?pe_lambda_grid). A user-supplied full grid also serves the reduced
# model unless a reduced grid was given, matching the reference
# implementation's `lamb_grid0 = lamb_grid`. Binary and count outcomes refit
# no reduced model and use only `full`.
#' @keywords internal
#' @noRd
.pe_default_grid <- function(n, p, outcome, model = "full") {
  if (outcome == "continuous") pe_lambda_grid(n, p, outcome, model)
  else seq(0.05, 1, length.out = 20)
}
.pe_resolve_grids <- function(lambda_grid, lambda_grid_reduced, n, p, outcome) {
  if (is.null(lambda_grid)) {
    full <- .pe_default_grid(n, p, outcome, "full")
    reduced <- if (is.null(lambda_grid_reduced))
      .pe_default_grid(n, p, outcome, "reduced") else lambda_grid_reduced
  } else {
    full <- lambda_grid
    reduced <- if (is.null(lambda_grid_reduced)) lambda_grid else lambda_grid_reduced
  }
  .validate_lambda_grid(full, "lambda_grid")
  .validate_lambda_grid(reduced, "lambda_grid_reduced")
  list(full = full, reduced = reduced)
}

# The target error rate for the mediator-screening step (FWER for Bonferroni,
# FDR for BH / BY) must be a single number strictly inside (0, 1): 0 would
# screen out every mediator and 1 would screen out none, and neither endpoint
# is a meaningful target rate.
#' @keywords internal
#' @noRd
.validate_error_level <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
    stop("`error_level` must be a single number in (0, 1).", call. = FALSE)
  invisible(TRUE)
}

# The confidence level for the total-indirect-effect interval must be a single
# number strictly inside (0, 1).
#' @keywords internal
#' @noRd
.validate_conf_level <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
    stop("`conf_level` must be a single number in (0, 1).", call. = FALSE)
  invisible(TRUE)
}

# Shared structural validation of (X, Y, M, Z): conforming row counts, a
# full-rank unpenalized design [X, Z] (the Wald inference inverts cross
# products of this block, so a collinear column makes the variance estimate
# singular), and no constant mediator columns (a zero-variance mediator cannot
# carry a signal and breaks standardization). Returns the coerced pieces plus
# n, q, and s = ncol(Z). Y-type checks (continuous vs 0/1 vs count) are left to
# the family worker, since only it knows the outcome model.
#' @keywords internal
#' @noRd
.pe_validate_inputs <- function(X, Y, M, Z, drop_constant = FALSE) {
  X <- .as_numeric_matrix(X, "X")
  M <- .as_numeric_matrix(M, "M")
  Y <- .as_numeric_vector(Y, "Y")
  n <- length(Y)
  # Column names, when the user supplied them (a named matrix or a
  # data.frame), travel to the output: the active set, the mediator table,
  # and the print footer are then read without a lookup against colnames(M),
  # and the per-exposure rows are named after the exposure columns. Captured
  # before any constant column is dropped, so they index the original `M`.
  mediator_names <- .pe_usable_names(M)
  exposure_names <- .pe_usable_names(X)

  if (nrow(X) != n)
    stop("`X` must have the same number of rows as the length of `Y`.",
         call. = FALSE)
  if (nrow(M) != n)
    stop("`M` must have the same number of rows as the length of `Y`.",
         call. = FALSE)

  if (!is.null(Z)) {
    Z <- .as_numeric_matrix(Z, "Z")
    if (nrow(Z) != n)
      stop("`Z` must have the same number of rows as the length of `Y`.",
           call. = FALSE)
  }

  # Every value must be finite: a missing, NaN, or infinite entry would
  # otherwise surface as an opaque error deep inside glmnet or the coordinate
  # descent, or silently poison the standardization. The method has no
  # missing-data mechanism; complete cases are the user's decision.
  for (nm in c("X", "Y", "M", "Z")) {
    v <- get(nm)
    if (!is.null(v) && !all(is.finite(v))) {
      bad <- sum(!is.finite(v))
      stop(sprintf(paste0("`%s` contains %d missing, NaN, or infinite value%s. ",
                          "POEM has no missing-data mechanism; supply complete ",
                          "cases (for example with stats::complete.cases())."),
                   nm, bad, if (bad == 1L) "" else "s"), call. = FALSE)
    }
  }

  # At least two candidate mediators: the penalized selection, the screen,
  # and the multiplicity adjustment are all defined over a set of mediators.
  if (ncol(M) < 2L)
    stop("`M` must have at least two columns (candidate mediators); a single ",
         "mediator calls for an ordinary mediation model.", call. = FALSE)

  # The exposures and confounders must vary: a constant column cannot be
  # standardized and would make the unpenalized design singular after scaling.
  for (nm in c("X", "Z")) {
    v <- get(nm)
    if (!is.null(v)) {
      const <- which(apply(v, 2L, function(col) length(unique(col)) < 2L))
      if (length(const))
        stop(sprintf("`%s` has constant column%s: %s. Remove %s; an intercept is fitted where the model needs one.",
                     nm, if (length(const) == 1L) "" else "s",
                     paste(const, collapse = ", "),
                     if (length(const) == 1L) "it" else "them"), call. = FALSE)
    }
  }

  # The unpenalized design [X, Z] must be full column rank. The mediation
  # inference inverts cross products formed from these columns; a rank
  # deficiency (a duplicated or collinear covariate) would make those
  # inversions fail or return garbage, so we stop early with a clear message.
  W <- if (is.null(Z)) X else cbind(X, Z)
  if (qr(W)$rank < ncol(W))
    stop("The design matrix formed by `X` and `Z` is singular (collinear ",
         "columns). Please remove redundant variables.", call. = FALSE)

  # A mediator with zero variance carries no information and cannot be
  # standardized (its standard deviation is zero). By default this is an error;
  # with drop_constant = TRUE the constant columns are dropped (with a warning)
  # and `kept` records the surviving columns' positions in the original `M`, so
  # the caller can map reported mediator indices back to the user's columns.
  m_var <- apply(M, 2L, stats::var)
  kept <- NULL
  if (any(m_var == 0)) {
    const_cols <- which(m_var == 0)
    if (!isTRUE(drop_constant))
      stop("The mediator matrix `M` contains constant (zero-variance) ",
           "columns: ", paste(const_cols, collapse = ", "),
           ". Remove them, or pass `drop_constant = TRUE` to drop them ",
           "automatically.", call. = FALSE)
    if (length(const_cols) >= ncol(M))
      stop("Every column of the mediator matrix `M` is constant ",
           "(zero variance); there is nothing to test.", call. = FALSE)
    warning(length(const_cols), " constant (zero-variance) mediator column",
            if (length(const_cols) > 1L) "s" else "", " dropped: ",
            paste(const_cols, collapse = ", "), ".", call. = FALSE)
    kept <- which(m_var != 0)
    M <- M[, kept, drop = FALSE]
  }

  list(X = X, Y = Y, M = M, Z = Z, n = n, q = ncol(X),
       s = if (is.null(Z)) 0L else ncol(Z), kept = kept,
       mediator_names = mediator_names, exposure_names = exposure_names)
}

# Column names usable as labels: present, non-empty, and unique. Anything
# else (no names, a blank, a duplicate) falls back to column positions.
#' @keywords internal
#' @noRd
.pe_usable_names <- function(x) {
  nm <- colnames(x)
  if (is.null(nm) || !all(nzchar(nm)) || anyDuplicated(nm)) return(NULL)
  nm
}

# Validation of a simulation design (simulate_mediation_data()). Each
# argument is checked by name so a bad value stops with a clear message
# instead of surfacing from the covariance construction or rnorm(). Not
# exported.
#' @keywords internal
#' @noRd
.validate_sim_design <- function(n, p, q, d, c1, c2, rho, sigma_y) {
  whole <- function(v, nm, at_least) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < at_least ||
        v != round(v))
      stop(sprintf("`%s` must be a single whole number of at least %d.",
                   nm, at_least), call. = FALSE)
  }
  scalar <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
      stop(sprintf("`%s` must be a single finite number.", nm), call. = FALSE)
  }
  whole(n, "n", 2L); whole(p, "p", 2L); whole(q, "q", 1L); whole(d, "d", 0L)
  scalar(c1, "c1"); scalar(c2, "c2"); scalar(rho, "rho"); scalar(sigma_y, "sigma_y")
  if (abs(rho) >= 1)
    stop("`rho` must lie strictly between -1 and 1.", call. = FALSE)
  if (sigma_y <= 0)
    stop("`sigma_y` must be positive.", call. = FALSE)
  invisible(TRUE)
}

# Standardize covariates and mediators (and center the continuous outcome)
# before fitting, following the preprocessing described in Yu and Kelley (in press): the
# penalty treats every mediator coefficient on a common scale, which requires
# the columns to share a scale. X, M, and Z are standardized to mean 0 and
# unit variance. For a continuous (Gaussian) outcome Y is centered so the
# intercept can be dropped from the linear model; for binary and count
# outcomes Y is a 0/1 or count response and is left untouched, since centering
# would destroy its meaning under the link.
#' @keywords internal
#' @noRd
.pe_scale_data <- function(X, Y, M, Z, center_y) {
  X <- base::scale(X)
  M <- base::scale(M)
  if (!is.null(Z)) Z <- base::scale(Z)
  if (center_y) Y <- Y - mean(Y)
  list(X = X, Y = Y, M = M, Z = Z)
}

# The global test rejects when the test statistic exceeds the upper-alpha
# quantile of a chi-square with q degrees of freedom (q = number of
# exposures), so the p-value is the chi-square upper-tail probability. Both the
# benchmark Wald statistic S_n and the power-enhanced statistic M_PE use this
# same reference distribution under the global null (Theorems on the limiting
# null distribution). lower.tail = FALSE is used rather than 1 - pchisq() so
# the far upper tail (p-values near machine zero, which the PE statistic
# routinely produces) keeps its precision instead of canceling to 0.
#' @keywords internal
#' @noRd
.pe_chisq_pvalue <- function(stat, df) {
  stats::pchisq(stat, df = df, lower.tail = FALSE)
}

# Optionally-parallel lapply used by the Monte Carlo and multi-group routines.
# On Unix with cores > 1 it forks via parallel::mclapply (parallel is a base
# package, so no extra dependency); otherwise it is an ordinary lapply. Callers
# that need reproducibility under forking set RNGkind("L'Ecuyer-CMRG") and a
# seed first, which gives parallel-safe streams. Not exported.
#' @keywords internal
#' @noRd
.pe_lapply <- function(X, FUN, cores = 1L) {
  if (cores > 1L && .Platform$OS.type == "unix")
    parallel::mclapply(X, FUN, mc.cores = cores)
  else
    lapply(X, FUN)
}
