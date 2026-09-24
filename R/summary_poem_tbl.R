#' Summarize a POEM table
#'
#' For a grouped table that compares the benchmark and power-enhanced tests
#' (any `poem_tbl` with a `group` column, a `pval_hdmm` column, and a
#' power-enhanced p-value column (`pval_pe` or `pval_pe_bonferroni`), such as
#' the output of [WHO_mediation_analysis()]), `summary()` digests it across
#' groups: how many groups each test detects mediation in, and in particular
#' the groups where the power-enhanced test detects mediation that the
#' benchmark misses, together with the most frequently flagged mediators. For
#' any other `poem_tbl` it falls back to the ordinary data-frame summary.
#'
#' @param object A `poem_tbl`.
#' @param alpha_level Significance level for counting a group as a detection.
#'   Default 0.05.
#' @param ... Ignored.
#'
#' @return For a grouped comparison table, an object of class
#'   `summary.poem_tbl` (a list, printed by its own method) with elements
#'   `outcome`, `alpha_level`, `n_groups`, `n_with_data`, `n_hdmm`, `n_pe` (groups
#'   detected by each test), `pe_only` (a `data.frame` of the groups the PE
#'   test detects but the benchmark does not), and `top_mediators` (a
#'   frequency table of flagged mediators). For any other `poem_tbl`, the
#'   value of [summary.data.frame()].
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [WHO_mediation_analysis()], [pe_mediation()].
#'
#' @family mediation tests
#'
#' @examples
#' # The global and income-group infant-mortality models on a 25-point
#' # tuning grid, which keeps the example fast; the article's tables use
#' # the default 100-point grid.
#' res <- WHO_mediation_analysis("imr", groupings = c("global", "income"),
#'                               lambda_grid = seq(0.1, 10, length.out = 25))
#' summary(res)
#'
#' @exportS3Method summary poem_tbl
summary.poem_tbl <- function(object, alpha_level = 0.05, ...) {
  .validate_alpha_level(alpha_level, "alpha_level")
  df <- as.data.frame(object)
  pe_col <- if ("pval_pe" %in% names(df)) "pval_pe"
            else if ("pval_pe_bonferroni" %in% names(df)) "pval_pe_bonferroni"
            else NA_character_
  # Only a grouped benchmark-vs-PE table gets the comparison digest; anything
  # else (a single fit, a power curve) falls back to the standard summary.
  if (!(all(c("group", "pval_hdmm") %in% names(df)) && !is.na(pe_col)))
    return(summary(df, ...))

  pe_p <- df[[pe_col]]
  pe_a <- if ("active_mediators" %in% names(df)) df$active_mediators
          else if ("active_bonferroni" %in% names(df)) df$active_bonferroni
          else rep(NA_character_, nrow(df))
  n_obs <- if ("n_countries" %in% names(df)) df$n_countries
           else rep(NA_integer_, nrow(df))
  has <- !is.na(df$pval_hdmm)                       # groups with data

  hdmm_sig <- has & df$pval_hdmm <= alpha_level
  pe_sig   <- has & pe_p <= alpha_level
  only     <- pe_sig & !hdmm_sig                    # PE detects, benchmark does not

  pe_only <- data.frame(
    group = df$group[only], n_countries = n_obs[only],
    pval_hdmm = df$pval_hdmm[only], pval_pe = pe_p[only],
    active_mediators = pe_a[only],
    row.names = NULL, stringsAsFactors = FALSE)

  codes <- unlist(strsplit(pe_a[pe_sig & pe_a != "none" & !is.na(pe_a)], ", "))
  top <- if (length(codes)) sort(table(codes), decreasing = TRUE) else integer(0)

  structure(list(
    outcome = attr(object, "outcome"),
    full_table = isTRUE(attr(object, "full_table")),
    alpha_level = alpha_level,
    n_groups = nrow(df), n_with_data = sum(has),
    n_hdmm = sum(hdmm_sig), n_pe = sum(pe_sig),
    pe_only = pe_only, top_mediators = top),
    class = "summary.poem_tbl")
}

#' @exportS3Method print summary.poem_tbl
print.summary.poem_tbl <- function(x, ...) {
  lab <- if (is.null(x$outcome)) "comparison" else toupper(x$outcome)
  cat(sprintf("POEM comparison: %s%s\n", lab,
              if (x$full_table) " (extended table)" else ""))
  cat(sprintf("  groups: %d (%d with data), alpha_level = %g\n",
              x$n_groups, x$n_with_data, x$alpha_level))
  cat(sprintf("  detected by benchmark (HDMM): %d   by power-enhanced (PE): %d\n",
              x$n_hdmm, x$n_pe))
  n_only <- nrow(x$pe_only)
  if (n_only == 0L) {
    cat("  PE detects no group the benchmark misses.\n")
  } else {
    cat(sprintf("  PE detects mediation in %d group%s the benchmark misses:\n",
                n_only, if (n_only > 1L) "s" else ""))
    show <- x$pe_only
    show$pval_hdmm <- signif(show$pval_hdmm, 3)
    show$pval_pe <- signif(show$pval_pe, 3)
    print(show, row.names = FALSE)
  }
  if (length(x$top_mediators)) {
    cat("  most-flagged mediators:\n")
    tm <- x$top_mediators
    for (i in seq_along(tm))
      cat(sprintf("    %-18s %d\n", names(tm)[i], tm[i]))
  }
  invisible(x)
}
