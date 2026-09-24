#' Power-enhanced mediation test for a continuous outcome
#'
#' Tests the global null hypothesis of no active mediator among a
#' high-dimensional set of candidate mediators, for a continuous
#' (linear-model) outcome, and reports both the benchmark Wald test on
#' the total indirect effect and the power-enhanced (PE) test that
#' remains powerful when individual mediation effects are heterogeneous
#' or contrasting. This is the continuous-outcome worker behind
#' [pe_mediation()]; call it directly when the outcome is continuous.
#'
#' @details
#' The model is the linear mediation pair
#' \deqn{Y = \alpha_m' M + \alpha_x' X + \alpha_z' Z + \varepsilon_y,
#'       \qquad M = \Gamma_x' X + \Gamma_z' Z + \varepsilon_m,}
#' with total indirect effect \eqn{\beta = \Gamma_x \alpha_m}. The
#' mediator coefficients \eqn{\alpha_m} are estimated by partial
#' penalized least squares with a SCAD penalty, the tuning parameter
#' chosen by the high-dimensional BIC. The benchmark statistic is the
#' Wald statistic \eqn{S_n = n \hat\beta' \hat\Sigma_\beta^{-1}
#' \hat\beta} of Guo et al. (2022); the power-enhanced statistic adds the
#' component \eqn{J_m} that accumulates the marginal signal from each
#' selected mediator (see [pe_mediation()] for the rationale and the
#' formula). Both are referred to a chi-square distribution with \eqn{q}
#' (the number of exposures) degrees of freedom under the global null.
#'
#' @inheritParams pe_mediation
#'
#' @return A tidy `data.frame` of class `poem_tbl` with rows
#'   `stat_hdmm` and `pval_hdmm` (the benchmark Wald test),
#'   `stat_pe`, `j_pe`, and `pval_pe` (the power-enhanced test and its
#'   PE component), `total_indirect_effect` with `total_indirect_lower`
#'   and `total_indirect_upper` (the estimate and a Wald interval at
#'   `conf_level`, the estimate plus or minus a standard normal quantile
#'   times its standard error; one row of each per exposure when `q > 1`,
#'   suffixed with the exposure's column name when `X` has column names
#'   and with `_1`, `_2`, ... otherwise), `n_active_mediators`, `df` (the
#'   chi-square degrees of freedom, equal to `q`), `n_candidate_mediators`,
#'   and `n_observations`. Attributes: `"active_mediators"` (the column
#'   positions in `M` identified as active; the print footer shows their
#'   column names when `M` has them), `"mediator_names"` and
#'   `"exposure_names"` (the column names of `M` and `X`, or `NULL` when
#'   unnamed), `"mediator_table"` (the per-mediator screening statistics,
#'   see [pe_mediators()]), `"method"`, `"error_level"`, `"conf_level"`,
#'   `"outcome"`, `"tuning"` (a list with the grid searched, the
#'   HBIC-selected `lambda_selected`, and `at_lower_end`, `TRUE` when the
#'   selection sat at the grid's smallest value; for a continuous outcome
#'   also the reduced model's grid and selection), `"empty_fit"`, and,
#'   with `report_all_methods = TRUE`, `"selection_by_method"` and
#'   `"pe_by_method"` (see [pe_selection()]).
#'
#'   When the penalized fit selects no mediator at any grid value there is
#'   nothing to test: the function warns (a condition of class
#'   `poem_empty_fit`), sets `"empty_fit"` to `TRUE`, and returns the same
#'   rows with both statistics 0, both p-values 1 (the article's Monte
#'   Carlo convention, under which an empty selection is a non-rejection),
#'   a total indirect effect of 0, and `NA` interval limits. A p-value of 1
#'   from such a fit is not evidence for the null; a grid reaching smaller
#'   values may select mediators.
#'
#' @references
#' Guo, X., Li, R., Liu, J., and Zeng, M. (2022). High-dimensional
#' mediation analysis for selecting DNA methylation loci mediating
#' childhood trauma and cortisol stress reactivity. \emph{Journal of the
#' American Statistical Association, 117}(539), 1110--1121.
#' \doi{10.1080/01621459.2022.2053136}
#'
#' Fan, J., and Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. \emph{Journal of the American
#' Statistical Association, 96}(456), 1348--1360.
#' \doi{10.1198/016214501753382273}
#'
#' Wang, L., Kim, Y., and Li, R. (2013). Calibrating nonconvex penalized
#' regression in ultra-high dimension. \emph{The Annals of Statistics,
#' 41}(5), 2505--2536. \doi{10.1214/13-AOS1159}
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()] for the front end, [pe_mediation_logistic()]
#'   and [pe_mediation_poisson()] for binary and count outcomes, and
#'   [simulate_mediation_data()] to generate data for trying the method.
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' # A contrasting setting: two active mediators whose indirect effects
#' # cancel, so the total indirect effect is zero. The Wald test is
#' # near powerless here; the PE test rejects and identifies active
#' # mediators.
#' d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
#'                              outcome = "continuous", c1 = 1)
#' pe_mediation_linear(d$X, d$Y, d$M)
#'
#' @export
pe_mediation_linear <- function(X, Y, M, Z = NULL,
                                method = c("Bonferroni", "BH", "BY"),
                                scale = TRUE, error_level = 0.05,
                                conf_level = 0.95, report_all_methods = FALSE,
                                drop_constant = FALSE,
                                lambda_grid = NULL,
                                lambda_grid_reduced = NULL) {
  method <- match.arg(method)
  .validate_error_level(error_level)
  .validate_conf_level(conf_level)
  # report_all_methods = TRUE screens with all three multiplicity methods (the
  # chosen one stays primary) so their active sets can be compared.
  methods <- if (isTRUE(report_all_methods))
    unique(c(method, "Bonferroni", "BH", "BY")) else method
  inp <- .pe_validate_inputs(X, Y, M, Z, drop_constant = drop_constant)
  X <- inp$X; Y <- inp$Y; M <- inp$M; Z <- inp$Z
  n <- inp$n; q <- inp$q; p <- ncol(M)

  # The tuning grids: the article's grids rescaled to this design unless the
  # user supplied one. A user-supplied full grid also serves the reduced model
  # unless a reduced grid was given, as in the reference implementation.
  grids <- .pe_resolve_grids(lambda_grid, lambda_grid_reduced, n, p, "continuous")
  lambda_grid <- grids$full; lambda_grid_reduced <- grids$reduced

  # A continuous outcome should not be binary; if it looks like counts, warn
  # but proceed (the user may genuinely want the linear model).
  if (length(unique(Y)) <= 2)
    stop("`Y` appears to be binary. Use pe_mediation_logistic() instead.",
         call. = FALSE)
  if (all(Y >= 0) && all(abs(Y - round(Y)) < .Machine$double.eps^0.5))
    warning("`Y` looks like a count variable; consider pe_mediation_poisson().",
            call. = FALSE)

  if (!is.logical(scale) || length(scale) != 1L || is.na(scale))
    stop("`scale` must be TRUE or FALSE.", call. = FALSE)
  if (scale) {
    sc <- .pe_scale_data(X, Y, M, Z, center_y = TRUE)
    X <- sc$X; Y <- sc$Y; M <- sc$M; Z <- sc$Z
  }
  S <- Z

  # --- Penalized fit + HBIC tuning for the full model. The grid of lambdas is
  # scored by HBIC and the minimizer (the last, if tied, matching the reference
  # implementation) is taken. The penalized design [M, X, S] does not change
  # across lambdas, so build it once and reuse it for every HBIC evaluation. ---
  MV_full <- if (length(S) == 0) cbind(M, X) else cbind(M, X, S)
  results <- lapply(lambda_grid, .HBIC_calc_linear, xx = X, yy = Y, mm = M,
                    S = S, n_imp = 0, MV = MV_full)
  hbic <- vapply(results, function(r) as.numeric(r$BIC), numeric(1))
  id <- utils::tail(which(hbic == min(hbic)), 1)
  res_full <- results[[id]]
  alpha0_hat <- res_full$alpha0; alpha1_hat <- res_full$alpha1
  alpha2_hat <- res_full$alpha2

  # --- Reduced-model fit (mediators on the outcome with only the intercept or
  # the confounders), used by the inference routine for the direct-effect
  # likelihood-ratio piece. Its design [M, intercept|S] is likewise constant. ---
  if (length(S) == 0) {
    intcpt <- matrix(rep(1, n), ncol = 1)
    MV_red <- cbind(M, intcpt)
    results0 <- lapply(lambda_grid_reduced, .HBIC_calc_linear, xx = intcpt,
                       yy = Y, mm = M, n_imp = 0, MV = MV_red)
  } else {
    MV_red <- cbind(M, S)
    results0 <- lapply(lambda_grid_reduced, .HBIC_calc_linear, xx = S,
                       yy = Y, mm = M, n_imp = 0, MV = MV_red)
  }
  hbic0 <- vapply(results0, function(r) as.numeric(r$BIC), numeric(1))
  id0 <- utils::tail(which(hbic0 == min(hbic0)), 1)
  alpha0_tld <- results0[[id0]]$alpha0
  alpha2_tld <- results0[[id0]]$alpha1
  tuning <- .pe_tuning_record(lambda_grid, lambda_grid[id],
                              lambda_grid_reduced, lambda_grid_reduced[id0])

  A     <- which(alpha0_hat != 0)    # selected mediators (support of alpha_m)
  A_tld <- which(alpha0_tld != 0)

  if (length(A) == 0L)
    return(.pe_relabel_mediators(
      .pe_empty_output(n, p, q, "continuous (linear)", method,
                       error_level, conf_level, methods = methods,
                       tuning = tuning,
                       exposure_names = inp$exposure_names),
      inp$kept, inp$mediator_names))

  M_A <- as.matrix(M[, A])
  if (length(A_tld) == 0L) {
    M_tld <- NULL; alpha0_tld_input <- rep(0, length(A))
  } else {
    M_tld <- as.matrix(M[, A_tld]); alpha0_tld_input <- alpha0_tld[A_tld]
  }

  # --- Benchmark Wald inference: S_n and the total indirect effect beta. ---
  hdmm <- .inference_linear(X, Y, M_A, S = S, M_tld = M_tld,
                            alpha0_hat = alpha0_hat[A],
                            alpha0_tld = alpha0_tld_input,
                            alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat,
                            alpha2_tld = alpha2_tld)
  stat_hdmm <- as.numeric(hdmm$Sn)
  pval_hdmm <- .pe_chisq_pvalue(stat_hdmm, df = q)

  # --- Confidence interval, power enhancement screen(s), and tidy output. ---
  .pe_relabel_mediators(
    .pe_finalize(X, Y, M_A, S, A, stat_hdmm, pval_hdmm, hdmm$beta_hat,
                 hdmm$var_beta, p, n, q, "gaussian", methods, error_level,
                 conf_level, "continuous (linear)", tuning = tuning,
                       exposure_names = inp$exposure_names),
    inp$kept, inp$mediator_names)
}
