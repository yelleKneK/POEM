#' Tidy printing for POEM result tables
#'
#' Every estimation and testing function in \pkg{POEM} returns a tidy
#' \code{data.frame} with a \code{term} column and one or more numeric
#' columns. A single numeric column routinely holds quantities on very
#' different scales: a whole-number count of mediators next to a
#' chi-square statistic in the hundreds next to a p-value of
#' \eqn{10^{-12}}. The base \code{\link[base]{print.data.frame}} method
#' formats a whole column with one common format, which forces either a
#' wall of trailing zeros or a slide into scientific notation. The
#' \code{poem_tbl} class supplies \code{print} and \code{format} methods
#' that format each value on its own terms: whole numbers (counts,
#' sample sizes, mediator indices) print without a decimal part, other
#' values print to a few significant figures, and p-values print to a
#' fixed number of decimal places with a \dQuote{< 0.0001} floor.
#'
#' The stored numeric values are never rounded. Only the display
#' changes, so any arithmetic you do on the returned object (a
#' confidence interval width, a further calculation) uses full
#' precision. To see more digits, raise \code{digits} or read the column
#' directly with \code{x$value}.
#'
#' This mirrors the display convention of the \pkg{DMAR} package, whose
#' house style \pkg{POEM} follows.
#'
#' @param x A \code{poem_tbl} object (a tidy \code{data.frame} returned
#'   by a POEM function).
#' @param digits Number of significant figures for non-integer values, a
#'   single positive whole number. Defaults to
#'   \code{getOption("poem.digits", 4L)}.
#' @param digits_p Number of decimal places for p-values, a single positive
#'   whole number. Defaults to 4. A p-value below \code{10^(-digits_p)}
#'   prints as \dQuote{< 0.0001}.
#' @param ... Additional arguments passed to
#'   \code{\link[base]{print.data.frame}}.
#'
#' @return \code{print.poem_tbl} returns \code{x} invisibly.
#'   \code{format.poem_tbl} returns a \code{data.frame} whose numeric
#'   columns have been formatted to character for display.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
#'                              outcome = "continuous")
#' fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#' fit                       # rounded for reading, with the p-value floor
#' fit$value[fit$term == "pval_pe"]   # full precision underneath
#'
#' @name poem_tbl
NULL

# Element-wise display formatter shared by format.poem_tbl / print.poem_tbl.
# A finite value that equals its own rounding (a whole number, such as a
# count of mediators or a sample size) prints with no decimal part; any
# other finite value prints to `digits` significant figures, letting base R
# choose fixed versus scientific notation per value (so a chi-square of 412.7
# reads plainly while a tiny variance reads as 1.2e-05). The whole-number
# branch is capped at 1e15 because beyond that doubles cannot represent
# integers exactly and "x == round(x)" stops being meaningful. Not exported.
.format_poem_value <- function(v, digits = 4L) {
  vapply(v, function(x) {
    if (is.na(x)) return("NA")
    if (is.finite(x) && x == round(x) && abs(x) < 1e15)
      return(format(x, scientific = FALSE, trim = TRUE))
    format(signif(x, digits), digits = digits, trim = TRUE)
  }, character(1L))
}

# Fixed-decimal p-value formatter. p-values read best at a fixed number of
# decimal places (the convention is four), never in scientific notation. A
# value below the smallest representable magnitude (10^(-digits_p)) prints as
# "< 0.0001" rather than rounding to "0.0000", which would misread as an exact
# zero. The PE tests can return p-values near machine zero (the chi-square
# upper tail at a statistic in the hundreds), so the floor matters here. Not
# exported.
.format_poem_pvalue <- function(p, digits_p = 4L) {
  thresh <- 10^(-digits_p)
  floor_label <- paste0("< ", formatC(thresh, format = "f", digits = digits_p))
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < thresh) return(floor_label)
    formatC(x, format = "f", digits = digits_p)
  }, character(1L))
}

# A digits argument must be a single positive whole number; anything else
# would make format() or formatC() error obscurely or print nothing.
.validate_digits <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 ||
      x != round(x))
    stop(sprintf("`%s` must be a single positive whole number.", name),
         call. = FALSE)
  invisible(TRUE)
}

# Tag a tidy result table as a poem_tbl so it prints via print.poem_tbl. Every
# tidy-returning POEM function routes its result through this helper just
# before returning it. `p_terms`, when supplied, names the rows of the long
# `value` column that hold p-values (so the format method gives them fixed
# decimals); a wide table with a literal `p_value` column needs no p_terms, as
# the format method finds that column by name. The poem_tbl class is inserted
# just before `data.frame`, after any leading tidy/glance subclass, so the
# subclass keeps dispatching its own methods while print falls through here.
# Idempotent. Not exported.
.as_poem_tbl <- function(out, p_terms = NULL) {
  if (!is.null(p_terms)) attr(out, "p_terms") <- p_terms
  cls <- class(out)
  if (!("poem_tbl" %in% cls))
    class(out) <- c(setdiff(cls, "data.frame"), "poem_tbl", "data.frame")
  out
}

#' @rdname poem_tbl
#' @export
format.poem_tbl <- function(x, digits = getOption("poem.digits", 4L),
                            digits_p = 4L, ...) {
  .validate_digits(digits, "digits")
  .validate_digits(digits_p, "digits_p")
  raw <- x
  class(raw) <- "data.frame"
  out <- raw
  num <- vapply(out, is.numeric, logical(1L))
  out[num] <- lapply(out[num], .format_poem_value, digits = digits)

  # The "primary" numeric column carries the per-term quantities in a long
  # table: `value`. The p_terms attribute names rows of this column whose
  # quantity is a p-value and so should print to fixed decimals.
  primary <- intersect("value", names(out)[num])
  primary <- if (length(primary)) primary[1L] else NA_character_

  # p-value columns in a wide table are detected by name: p_value, p.value,
  # pval_* (the fit and WHO tables), and screen_p (the mediator table).
  pcols <- grep("^(p_value|p\\.value|pval(_.*)?|screen_p)$", names(out)[num],
                value = TRUE)
  for (nm in pcols)
    out[[nm]] <- .format_poem_pvalue(raw[[nm]], digits_p = digits_p)

  # Named rows of the long-format `value` column are detected via p_terms.
  p_terms <- attr(x, "p_terms")
  if (!is.null(p_terms) && "term" %in% names(out) && !is.na(primary)) {
    idx <- raw$term %in% p_terms
    out[[primary]][idx] <- .format_poem_pvalue(raw[[primary]][idx],
                                               digits_p = digits_p)
  }
  out
}

#' @rdname poem_tbl
#' @export
print.poem_tbl <- function(x, digits = getOption("poem.digits", 4L),
                           digits_p = 4L, ...) {
  disp <- format(x, digits = digits, digits_p = digits_p)
  # row.names / right travel through `...` so a caller can override them, but
  # default to the tidy look (no row numbers, left-aligned).
  dots <- list(...)
  if (is.null(dots[["row.names"]])) dots[["row.names"]] <- FALSE
  if (is.null(dots[["right"]]))     dots[["right"]]     <- FALSE
  do.call(print.data.frame, c(list(disp), dots))

  # A short, human-readable footer naming the outcome model and any active
  # mediators the test identified. The numbers are already in the table; the
  # footer just makes the headline reading effortless.
  # A fit names its outcome model; a study or WHO table names its outcome.
  outcome <- attr(x, "outcome")
  if (!is.null(outcome))
    cat(sprintf("\n%s: %s\n", if ("term" %in% names(x)) "Outcome model" else "Outcome",
                outcome))
  active <- attr(x, "active_mediators")
  if (!is.null(active)) {
    if (length(active) == 0L) {
      cat("Active mediators identified: none\n")
    } else {
      # By column name when `M` had names, by column position otherwise.
      nm <- attr(x, "mediator_names")
      lab <- if (is.null(nm)) active else nm[active]
      cat(sprintf("Active mediators identified (%d): %s\n",
                  length(active), paste(lab, collapse = ", ")))
    }
  }
  # The tuning record: which lambda HBIC chose from which grid, and whether it
  # sat at the grid's lower end (where the floor, not HBIC, decided the fit).
  tun <- attr(x, "tuning")
  if (!is.null(tun)) {
    g <- tun$lambda_grid
    cat(sprintf("Tuning parameter (HBIC): lambda = %s from %d values in [%s, %s]%s\n",
                format(signif(tun$lambda_selected, 3)), length(g),
                format(signif(min(g), 3)), format(signif(max(g), 3)),
                if (isTRUE(tun$at_lower_end)) " (the grid's lower end)" else ""))
  }
  if (isTRUE(attr(x, "empty_fit")))
    cat("Penalized selection: empty at every lambda; no test was computed.\n")
  invisible(x)
}
