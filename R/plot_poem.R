#' Plot a POEM result
#'
#' A base-graphics plot for the two kinds of table POEM returns. For a
#' power curve from [pe_power_curve()] it draws the empirical rejection
#' rate of the benchmark and power-enhanced tests against the
#' signal-strength grid, with the nominal level marked. For a single
#' [pe_mediation()] fit it draws the per-mediator screening evidence
#' (\eqn{-\log_{10}} of the screening p-value) for each selected
#' candidate mediator, with the active ones highlighted, so it is clear
#' which mediators drove the result.
#'
#' @param x A `poem_tbl` from [pe_power_curve()] or [pe_mediation()].
#' @param ... Further graphical parameters passed to the underlying plot.
#'
#' @return `x`, invisibly. Called for the plot it draws.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_power_curve()], [pe_mediation()], [pe_mediators()].
#'
#' @family mediation simulation
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
#'                              pattern = "contrasting", c1 = 1)
#' fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#' plot(fit)
#'
#' @method plot poem_tbl
#' @export
plot.poem_tbl <- function(x, ...) {
  if (all(c("c1", "rejection_hdmm", "rejection_pe") %in% names(x))) {
    # Power-curve view: two rejection-rate curves against the signal grid.
    alpha_level <- attr(x, "alpha_level")
    plot(x$c1, x$rejection_pe, type = "b", pch = 19, ylim = c(0, 1),
         xlab = "signal strength (c1)", ylab = "empirical rejection rate",
         ...)
    graphics::lines(x$c1, x$rejection_hdmm, type = "b", pch = 1, lty = 2)
    if (!is.null(alpha_level)) graphics::abline(h = alpha_level, col = "grey60", lty = 3)
    graphics::legend("right", c("power-enhanced", "benchmark Wald"),
                     pch = c(19, 1), lty = c(1, 2), bty = "n")
    return(invisible(x))
  }

  tab <- attr(x, "mediator_table")
  if (!is.null(tab) && nrow(tab) > 0L) {
    # Per-mediator view: screening evidence, active mediators highlighted.
    ev <- -log10(pmax(tab$screen_p, .Machine$double.xmin))
    cols <- ifelse(tab$selected, "black", "grey70")
    plot(seq_len(nrow(tab)), ev, pch = 19, col = cols, xaxt = "n",
         xlab = "selected candidate mediator",
         ylab = expression(-log[10](screening~p)), ...)
    graphics::axis(1, at = seq_len(nrow(tab)),
                   labels = if ("name" %in% names(tab)) tab$name else tab$mediator)
    graphics::legend("topright", c("active", "screened out"),
                     pch = 19, col = c("black", "grey70"), bty = "n")
    return(invisible(x))
  }

  message("No POEM plot is defined for this table.")
  invisible(x)
}

#' Broom verbs for POEM results
#'
#' `tidy()` and `glance()` methods so a POEM result composes with the
#' broom ecosystem; the generics are re-exported, so `tidy(fit)` and
#' `glance(fit)` work after `library(POEM)` alone. For a [pe_mediation()]
#' fit, `tidy()` returns the per-mediator table (see [pe_mediators()]) and
#' `glance()` returns a one-row summary of the two tests and the total
#' indirect effect. For a [pe_power_curve()] table both return the table
#' itself (it is already one row per grid point).
#'
#' @param x A `poem_tbl` result.
#' @param ... Unused, for generic compatibility.
#'
#' @return A `data.frame`. For a fit, `tidy()` has the columns of
#'   [pe_mediators()] and `glance()` the columns `stat_hdmm`, `pval_hdmm`,
#'   `stat_pe`, `pval_pe`, `total_indirect_effect` (one column per
#'   exposure when there are several, suffixed with the exposure's column
#'   name or, for an unnamed `X`, `_1`, `_2`, ...), `n_active_mediators`,
#'   and `n_observations`.
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 120, p = 40, pattern = "contrasting",
#'                              outcome = "continuous", c1 = 1)
#' fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#' tidy(fit)
#' glance(fit)
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()], [pe_mediators()].
#'
#' @name poem_broom
#' @importFrom generics tidy
#' @export
tidy.poem_tbl <- function(x, ...) {
  tab <- attr(x, "mediator_table")
  if (!is.null(tab)) return(as.data.frame(tab))
  as.data.frame(unclass(x))
}

#' @rdname poem_broom
#' @importFrom generics glance
#' @export
glance.poem_tbl <- function(x, ...) {
  if (!"term" %in% names(x)) return(as.data.frame(unclass(x)))
  v <- function(t) {
    hit <- x$value[x$term == t]
    if (length(hit)) hit[1] else NA_real_
  }
  # One total-indirect-effect column per exposure (the rows are suffixed
  # _1, _2, ... when there are several exposures).
  beta_terms <- grep("^total_indirect_effect(_.+)?$", x$term, value = TRUE)
  beta <- stats::setNames(as.list(x$value[match(beta_terms, x$term)]), beta_terms)
  cbind(data.frame(stat_hdmm = v("stat_hdmm"), pval_hdmm = v("pval_hdmm"),
                   stat_pe = v("stat_pe"), pval_pe = v("pval_pe"),
                   stringsAsFactors = FALSE),
        as.data.frame(beta, stringsAsFactors = FALSE),
        data.frame(n_active_mediators = v("n_active_mediators"),
                   n_observations = v("n_observations"),
                   stringsAsFactors = FALSE))
}

#' @rdname poem_broom
#' @export
generics::tidy

#' @rdname poem_broom
#' @export
generics::glance
