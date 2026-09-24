# Assemble the tidy result table shared by all three pe_mediation_*() workers,
# so the continuous, binary, and count outcomes return one consistent schema.
# The table reports, side by side, the benchmark Wald test on the total
# indirect effect (labelled "hdmm", after Guo et al.) and the power-enhanced
# test ("pe"), as the article's data-analysis tables
# do, plus a confidence interval for the total indirect effect. The active
# mediator set, the per-mediator table, and any per-method selections are
# variable length and not numeric, so they travel as attributes (keeping the
# value column numeric) and are surfaced through the print footer and the
# accessors pe_mediators() / pe_selection(). Not exported.
#
# Arguments are the already-computed scalars from a worker:
#   stat_hdmm, pval_hdmm   the Wald statistic S_n and its chi-square_q p-value
#   stat_pe, j_pe, pval_pe the PE statistic M_PE = S_n + J_m, J_m, and p-value
#   beta_hat               the estimated total indirect effect (length q)
#   beta_lower, beta_upper the confidence limits for beta (length q), or NULL
#   selected_mediators     indices (into 1..p) of the identified active mediators
#   n, p, q                sample size, number of candidate mediators, exposures
#   outcome_label          human-readable outcome model name for the footer
#   method, error_level    the multiplicity method and its target FWER/FDR level
#   conf_level             the confidence level for the beta interval
#   mediator_table         per-mediator data.frame (attribute), or NULL
#   selection_by_method    named list of active sets per method (attribute), or NULL
#   pe_by_method           data.frame of the PE statistic/p-value per method
#                          (attribute), or NULL; columns method, stat_pe, j_pe,
#                          pval_pe, n_active
#' @keywords internal
#' @noRd
.pe_assemble_output <- function(stat_hdmm, pval_hdmm, stat_pe, j_pe, pval_pe,
                                beta_hat, beta_lower = NULL, beta_upper = NULL,
                                selected_mediators, n, p, q,
                                outcome_label, method, error_level,
                                conf_level = 0.95, mediator_table = NULL,
                                selection_by_method = NULL,
                                pe_by_method = NULL, tuning = NULL,
                                empty_fit = FALSE, exposure_names = NULL) {
  beta_hat <- as.numeric(beta_hat)
  # One total-indirect-effect row when there is a single exposure (the common
  # case), otherwise one row per exposure so the value column stays scalar,
  # suffixed with the exposure's column name when `X` has usable names and
  # with its position otherwise.
  if (length(beta_hat) == 1L) {
    beta_terms <- "total_indirect_effect"
    lower_terms <- "total_indirect_lower"
    upper_terms <- "total_indirect_upper"
  } else {
    sfx <- if (!is.null(exposure_names) &&
               length(exposure_names) == length(beta_hat))
      exposure_names else seq_along(beta_hat)
    beta_terms  <- paste0("total_indirect_effect_", sfx)
    lower_terms <- paste0("total_indirect_lower_", sfx)
    upper_terms <- paste0("total_indirect_upper_", sfx)
  }

  # The confidence-interval rows are always present so the row schema never
  # changes; the empty-fit path supplies NA limits (no interval exists when no
  # mediator was selected).
  if (is.null(beta_lower)) beta_lower <- rep(NA_real_, length(beta_hat))
  if (is.null(beta_upper)) beta_upper <- rep(NA_real_, length(beta_hat))
  terms  <- c("stat_hdmm", "pval_hdmm", "stat_pe", "j_pe", "pval_pe",
              beta_terms, lower_terms, upper_terms,
              "n_active_mediators", "df", "n_candidate_mediators",
              "n_observations")
  values <- c(stat_hdmm, pval_hdmm, stat_pe, j_pe, pval_pe, beta_hat,
              as.numeric(beta_lower), as.numeric(beta_upper),
              length(selected_mediators), q, p, n)

  out <- data.frame(term = terms, value = values, stringsAsFactors = FALSE)

  # p-value rows print to fixed decimals via the poem_tbl format method.
  out <- .as_poem_tbl(out, p_terms = c("pval_hdmm", "pval_pe"))
  attr(out, "outcome")             <- outcome_label
  attr(out, "active_mediators")    <- selected_mediators
  attr(out, "method")              <- method
  attr(out, "error_level")         <- error_level
  attr(out, "conf_level")          <- conf_level
  attr(out, "mediator_table")      <- mediator_table
  attr(out, "selection_by_method") <- selection_by_method
  attr(out, "pe_by_method")        <- pe_by_method
  # The tuning record: the grid(s) searched and the HBIC-selected value(s),
  # so a reader can see which lambda produced the fit and whether it sat at
  # the grid's lower end (see pe_lambda_grid()).
  attr(out, "tuning")              <- tuning
  attr(out, "empty_fit")           <- isTRUE(empty_fit)
  attr(out, "exposure_names")      <- exposure_names
  out
}

# The tuning record attached to every fit. `grid` and `selected` are the
# candidate grid and the HBIC-selected value for the full model; the reduced
# model's pair is recorded when the family refits one (continuous outcome).
# `at_lower_end` flags a selection at the smallest grid value, where the grid's
# floor rather than HBIC determined the fit. Not exported.
#' @keywords internal
#' @noRd
.pe_tuning_record <- function(grid, selected, grid_reduced = NULL,
                              selected_reduced = NULL) {
  rec <- list(lambda_grid = grid, lambda_selected = selected,
              at_lower_end = isTRUE(selected <= min(grid)),
              at_upper_end = isTRUE(selected >= max(grid)))
  if (!is.null(grid_reduced)) {
    rec$lambda_grid_reduced <- grid_reduced
    rec$lambda_selected_reduced <- selected_reduced
  }
  rec
}

# The warning raised when the penalized fit selects no mediator at any grid
# value. It carries the class "poem_empty_fit" so the Monte Carlo functions
# can count these fits without printing hundreds of warnings. Not exported.
#' @keywords internal
#' @noRd
.pe_empty_fit_warning <- function(grid) {
  msg <- sprintf(paste0(
    "No mediator survived the penalized selection at any of the %d values of ",
    "`lambda_grid` (%.3g to %.3g); the benchmark and power-enhanced tests were ",
    "not computed and are reported as 0 with p-values of 1, the convention of ",
    "Yu and Kelley (in press) for an empty selection. A grid extending to ",
    "smaller values may select mediators; see ?pe_lambda_grid."),
    length(grid), min(grid), max(grid))
  warning(structure(class = c("poem_empty_fit", "warning", "condition"),
                    list(message = msg, call = NULL)))
}

# Post-inference finalization shared by the three workers. Given the benchmark
# Wald pieces (statistic, p-value, beta, and its asymptotic covariance) and the
# selected mediators, it (1) forms a Wald confidence interval for the total
# indirect effect, (2) runs the power-enhancement screen for each requested
# multiplicity method (the first is primary; the rest, if any, are recorded for
# comparison), and (3) assembles the tidy output. Centralizing this keeps the
# workers identical past the family-specific fit. Not exported.
#
# The covariance var_beta is the asymptotic covariance of sqrt(n)*beta_hat (the
# matrix the Wald statistic S_n = n beta' var_beta^{-1} beta inverts), so the
# variance of beta_hat itself is var_beta / n and the standard error is the
# square root of its diagonal over n.
#' @keywords internal
#' @noRd
.pe_finalize <- function(X, Y, M_A, S, A, stat_hdmm, pval_hdmm, beta_hat,
                         var_beta, p, n, q, y_family, methods, error_level,
                         conf_level, outcome_label, tuning = NULL,
                         exposure_names = NULL) {
  vb <- as.matrix(var_beta)
  se <- sqrt(pmax(diag(vb), 0) / n)
  z  <- stats::qnorm(1 - (1 - conf_level) / 2)
  beta_lower <- as.numeric(beta_hat) - z * se
  beta_upper <- as.numeric(beta_hat) + z * se

  # Screen with each requested method off the SAME penalized fit (the penalized
  # estimation, the expensive part, is already done; the screen only refits the
  # cheap path regressions). The first method drives the primary output.
  pe_list <- lapply(methods, function(m)
    .pe_component(X, Y, M_A, S, A, Sn = stat_hdmm, p_total = p, n = n, q = q,
                  y_family = y_family, method = m, error_level = error_level))
  names(pe_list) <- methods
  primary <- pe_list[[1L]]
  selection_by_method <- if (length(methods) > 1L)
    lapply(pe_list, function(pc) pc$selected_mediators) else NULL

  # The PE statistic, J_m, and p-value differ by multiplicity method (each
  # method screens a different active set into J_m), so record all of them.
  # This is what pe_selection() reports side by side and the article's
  # "full" data-analysis tables show as the PE / PE_BH / PE_BY columns.
  pe_by_method <- data.frame(
    method   = methods,
    stat_pe  = vapply(pe_list, function(pc) pc$stat_pe, numeric(1)),
    j_pe     = vapply(pe_list, function(pc) pc$j_pe, numeric(1)),
    pval_pe  = vapply(pe_list, function(pc) pc$pval_pe, numeric(1)),
    n_active = vapply(pe_list, function(pc) length(pc$selected_mediators),
                      integer(1)),
    row.names = NULL, stringsAsFactors = FALSE)

  .pe_assemble_output(stat_hdmm, pval_hdmm, primary$stat_pe, primary$j_pe,
                      primary$pval_pe, beta_hat, beta_lower, beta_upper,
                      primary$selected_mediators, n, p, q, outcome_label,
                      methods[1L], error_level, conf_level,
                      mediator_table = primary$mediator_table,
                      selection_by_method = selection_by_method,
                      pe_by_method = pe_by_method, tuning = tuning,
                      exposure_names = exposure_names)
}

# Map mediator indices in a finished result from positions in the (possibly
# reduced) mediator matrix back to positions in the user's original `M`, and
# attach the mediator names when `M` had them. `kept` is the vector of
# surviving original-column indices when drop_constant dropped zero-variance
# columns (a reduced index j becomes kept[j]); NULL when nothing was dropped.
# `mediator_names` are the column names of the user's original `M` (NULL when
# unnamed); with names, the mediator table gains a `name` column and the
# names travel in the "mediator_names" attribute for the print footer,
# pe_selection(), and plot(). Not exported.
#' @keywords internal
#' @noRd
.pe_relabel_mediators <- function(out, kept, mediator_names = NULL) {
  if (!is.null(kept)) {
    am <- attr(out, "active_mediators")
    if (length(am)) attr(out, "active_mediators") <- kept[am]
    sbm <- attr(out, "selection_by_method")
    if (!is.null(sbm))
      attr(out, "selection_by_method") <-
        lapply(sbm, function(idx) if (length(idx)) kept[idx] else idx)
    mt <- attr(out, "mediator_table")
    if (!is.null(mt) && nrow(mt)) {
      mt$mediator <- kept[mt$mediator]
      attr(out, "mediator_table") <- mt
    }
  }
  attr(out, "mediator_names") <- mediator_names
  mt <- attr(out, "mediator_table")
  if (!is.null(mediator_names) && !is.null(mt)) {
    rest <- setdiff(names(mt), "mediator")
    mt <- cbind(mt[, "mediator", drop = FALSE],
                data.frame(name = mediator_names[mt$mediator],
                           stringsAsFactors = FALSE),
                mt[, rest, drop = FALSE])
    rownames(mt) <- NULL
    attr(out, "mediator_table") <- mt
  }
  out
}

# The result returned when the penalized fit selects no candidate-active
# mediator at any grid value (the selected set S-hat is empty). With no
# selected mediator there is no indirect effect to test: J_m = 0 by definition,
# and the benchmark statistic is reported as 0 with a p-value of 1, the
# convention of the article's Monte Carlo accounting (an empty selection is a
# non-rejection). The result keeps the full row schema (the interval rows are
# NA), carries `empty_fit = TRUE`, and is announced by a classed warning so it
# cannot be mistaken for evidence for the null.
#
# `methods` is the vector of multiplicity methods requested (length > 1 when
# report_all_methods = TRUE); with no selected mediators every method gives the
# same all-zero result, so they share J_m = 0, stat_pe = 0, pval_pe = 1.
#' @keywords internal
#' @noRd
.pe_empty_output <- function(n, p, q, outcome_label, method, error_level,
                             conf_level = 0.95, methods = method,
                             tuning = NULL, exposure_names = NULL) {
  pe_by_method <- data.frame(
    method = methods, stat_pe = 0, j_pe = 0, pval_pe = 1, n_active = 0L,
    row.names = NULL, stringsAsFactors = FALSE)
  selection_by_method <- if (length(methods) > 1L)
    stats::setNames(rep(list(integer(0)), length(methods)), methods) else NULL
  if (!is.null(tuning)) .pe_empty_fit_warning(tuning$lambda_grid)
  .pe_assemble_output(
    stat_hdmm = 0, pval_hdmm = 1, stat_pe = 0, j_pe = 0, pval_pe = 1,
    beta_hat = rep(0, q), beta_lower = NULL, beta_upper = NULL,
    selected_mediators = integer(0), n = n, p = p, q = q,
    outcome_label = outcome_label, method = methods[1L],
    error_level = error_level, conf_level = conf_level,
    mediator_table = data.frame(
      mediator = integer(0), t_outcome = numeric(0), t_exposure = numeric(0),
      screen_p = numeric(0), selected = logical(0)),
    selection_by_method = selection_by_method,
    pe_by_method = pe_by_method, tuning = tuning, empty_fit = TRUE,
    exposure_names = exposure_names)
}
