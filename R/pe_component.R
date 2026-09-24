# The power enhancement (PE) component and active-mediator selection.
#
# This is the mechanism that distinguishes the PE tests from a plain Wald test
# on the total indirect effect, and it is identical across the continuous,
# binary, and count outcome models (only the outcome regression that supplies
# the mediator-on-outcome statistics changes family). The prior PEmediation
# package inlined a byte-for-byte copy of this block in each of its three
# workers; here it lives once.
#
# The idea (Fan, Liao, and Yao, 2015, adapted to mediation). After the
# penalized fit selects a small set A of candidate-active mediators, we form,
# for each selected mediator j and each exposure i, two t-like statistics:
#
#   t_alpha[j]   from the outcome regression  Y ~ M_A + X (+ Z): the M -> Y path
#   t_Gamma[i,j] from the mediator regression M_A ~ X (+ Z):     the X -> M path
#
# A mediator is "active" only if BOTH paths are nonzero, so the individual
# mediation signal is screened on the LARGER of the two p-values,
# max(p_alpha, p_Gamma): a mediator passes only when both paths are
# significant. The PE component accumulates the product of magnitudes
# |t_alpha * t_Gamma| over the pairs that pass the screen, scaled by sqrt(p):
#
#   J_m = sqrt(p) * sum over passing (i,j) of |t_alpha[j] * t_Gamma[i,j]|.
#
# Because J_m sums magnitudes, two mediators whose indirect effects have
# opposite signs REINFORCE in J_m instead of canceling, which is exactly the
# heterogeneous / contrasting case where the total-indirect-effect Wald test
# S_n has no power. The PE test statistic is M_PE = S_n + J_m. Under the global
# null of no active mediator the screen rejects everything with probability
# approaching one, so J_m = 0 and M_PE has the same chi-square_q limit as S_n;
# the size of the test is preserved while power is enhanced.
#
# The screening threshold uses a Bonferroni correction over the s * q tested
# pairs (s = |A|), divided by an extra log(log n). That slowly diverging factor
# drives the probability of selecting a spurious mediator under the null to
# zero (it can be any diverging sequence; log log n is among the slowest, which
# avoids needless conservativeness). The Benjamini-Hochberg and
# Benjamini-Yekutieli alternatives trade the familywise error rate guarantee
# for false discovery rate control, which is less conservative when many
# mediators are active.
#
# Inputs:
#   X, Y, M_A   the (scaled) exposures, outcome, and SELECTED mediator columns
#   S           the (scaled) confounders Z, or NULL
#   A           integer indices (into the full 1..p) of the selected mediators
#   Sn          the benchmark Wald statistic from the family inference routine
#   p_total     p, the number of CANDIDATE mediators (the sqrt(p) scale factor)
#   n, q        sample size and number of exposures
#   y_family    "gaussian", "binomial", or "poisson"; the outcome glm family
#   method      "Bonferroni" (FWER), "BH", or "BY" (FDR)
#   error_level the target FWER (Bonferroni) or FDR (BH / BY) level
#
# Returns a list: stat_pe (M_PE), j_pe (J_m), pval_pe, and selected_mediators
# (the active set, as indices into 1..p).
#' @keywords internal
#' @noRd
.pe_component <- function(X, Y, M_A, S, A, Sn, p_total, n, q,
                          y_family, method, error_level) {
  k <- length(A)   # number of penalized-selected (candidate-active) mediators

  # ------------------------------------------------------------------
  # Path M -> Y: t-statistics for the selected mediators' outcome
  # coefficients, from the outcome regression on the selected mediators,
  # exposures, and confounders. For a Gaussian outcome glm() reproduces lm()
  # and column 3 is the t value; for binomial / Poisson it is the z value.
  # Either way it is the standardized coefficient alpha_m,j / se. The data
  # frame is built so that "Y ~ ." expands to the mediator columns first,
  # hence the M coefficients occupy rows 2:(k+1) (row 1 is the intercept).
  # ------------------------------------------------------------------
  if (is.null(S)) {
    ymx <- data.frame(Y = Y, M = M_A, X = X)
  } else {
    ymx <- data.frame(Y = Y, M = M_A, X = X, S = S)
  }
  res <- summary(stats::glm(Y ~ ., family = y_family, data = ymx))$coefficients
  t_alpha <- res[2:(k + 1), 3]                 # length k
  p_alpha <- 2 * (1 - stats::pnorm(abs(t_alpha)))

  # ------------------------------------------------------------------
  # Path X -> M: t-statistics for the exposure coefficients in the mediator
  # regression M_A ~ X (+ Z), one per (exposure i, mediator j) pair. With more
  # than one selected mediator the response is multivariate, so summary()
  # returns one summary.lm per mediator column and we pull the X rows (rows
  # 2:(1+q), after the intercept) from each. The matrices are referenced
  # directly (the regression is multivariate least squares of the matrix M_A on
  # the matrix X), which is what the formula resolves to.
  # ------------------------------------------------------------------
  if (is.null(S)) {
    res_mx <- summary(stats::lm(M_A ~ X))
  } else {
    res_mx <- summary(stats::lm(M_A ~ X + S))
  }
  if (k == 1L) {
    est_Gamma   <- res_mx$coefficients[2:(1 + q), 1]
    t_Gamma     <- res_mx$coefficients[2:(1 + q), 3]
  } else {
    t_Gamma     <- sapply(seq_len(k),
                          function(i) res_mx[[i]]$coefficients[2:(1 + q), 3])
  }
  p_Gamma <- 2 * (1 - stats::pnorm(abs(t_Gamma)))

  # Shape every X -> M quantity as a q-by-k matrix (rows = exposures, columns =
  # selected mediators). When q == 1 the sapply / single-fit results come back
  # as length-k vectors, so transpose into a 1-by-k row.
  if (q == 1L) {
    t_Gamma <- as.matrix(t(t_Gamma))
    p_Gamma <- as.matrix(t(p_Gamma))
  } else {
    t_Gamma <- as.matrix(t_Gamma)
    p_Gamma <- as.matrix(p_Gamma)
  }

  # The M -> Y statistics depend on j only, so tile them across the q exposure
  # rows to align with the q-by-k X -> M matrices for the elementwise product.
  t_alpha_mat <- matrix(t_alpha, nrow = q, ncol = k, byrow = TRUE)
  p_alpha_mat <- matrix(p_alpha, nrow = q, ncol = k, byrow = TRUE)

  # The per-pair screening p-value is the LARGER of the two path p-values: a
  # pair survives only when both the X -> M and the M -> Y path are significant.
  screen_p <- pmax(p_alpha_mat, p_Gamma)

  # ------------------------------------------------------------------
  # Multiplicity screen, then J_m and the active-mediator set. The three
  # methods differ only in how the per-pair threshold is formed; the J_m sum
  # and the selection rule are the same. selected_mediators collapses the
  # q-by-k indicator to the mediators (columns) that pass for at least one
  # exposure (row), mapped back to the original 1..p indices via A.
  # ------------------------------------------------------------------
  if (method == "Bonferroni") {
    # Familywise error rate control by Bonferroni over the k * q pairs, with
    # the extra log(log n) that vanishes the spurious-selection probability.
    threshold <- error_level / (k * q * log(log(n)))
    indicator <- (screen_p <= threshold)
  } else {
    # BH or BY: convert the screening p-values to adjusted p-values (q-values)
    # over the k * q pairs, then compare to error_level / log(log n). BH is
    # valid under independence or positive dependence among mediators; BY is
    # valid under arbitrary dependence (at the cost of more conservatism).
    adj <- stats::p.adjust(screen_p, method = method, n = k * q)
    indicator <- (adj <= (error_level / log(log(n))))
    indicator <- matrix(indicator, nrow = q, ncol = k)
  }

  # J_m: sqrt(p) times the sum of |t_alpha * t_Gamma| over the surviving pairs.
  j_pe <- sqrt(p_total) * sum(abs((t_alpha_mat * t_Gamma)[which(indicator)]))

  # Active mediators: those columns (mediators) flagged for some exposure.
  if (q == 1L) {
    selected_logical <- as.logical(indicator)
  } else {
    selected_logical <- apply(indicator, 2L, any)
  }
  selected_mediators <- A[selected_logical]

  # Per-mediator summary, one row per penalized-selected mediator (the set A):
  # its standardized M -> Y statistic, its strongest standardized X -> M
  # statistic across exposures, the smallest screening p-value across
  # exposures, and whether it survived the screen. This is what pe_mediators()
  # surfaces so a user can see why each mediator was or was not flagged.
  mediator_table <- data.frame(
    mediator = A,
    t_outcome  = t_alpha,                                  # M -> Y, per mediator
    t_exposure = if (q == 1L) as.numeric(t_Gamma)
                 else apply(abs(t_Gamma), 2L, max),        # strongest X -> M
    screen_p   = if (q == 1L) as.numeric(screen_p)
                 else apply(screen_p, 2L, min),
    selected   = selected_logical,
    stringsAsFactors = FALSE)

  stat_pe <- Sn + j_pe
  list(stat_pe = stat_pe,
       j_pe = j_pe,
       pval_pe = .pe_chisq_pvalue(stat_pe, df = q),
       selected_mediators = selected_mediators,
       mediator_table = mediator_table)
}
