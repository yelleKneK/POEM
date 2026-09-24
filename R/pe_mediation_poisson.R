#' Power-enhanced mediation test for a count outcome
#'
#' Tests the global null hypothesis of no active mediator for a count
#' outcome modeled with Poisson (log-link) regression, reporting both
#' the benchmark Wald test on the total indirect effect and the
#' power-enhanced (PE) test. This is the count-outcome worker behind
#' [pe_mediation()].
#'
#' @details
#' The outcome follows a Poisson mediation model,
#' \eqn{Y \mid M, X, Z \sim \mathrm{Poisson}(\exp(\alpha_m' M +
#' \alpha_x' X + \alpha_z' Z))}, with a linear mediator model as in
#' [pe_mediation_linear()]. Estimation, tuning, and the power
#' enhancement construction follow the generalized-outcome framework of
#' Guo et al. (2024); the penalized Poisson fit uses a glmnet lasso warm
#' start followed by a weighted local linear approximation step. The
#' total indirect effect \eqn{\beta = \Gamma_x \alpha_m} is on the
#' log-mean (link) scale, not the count scale.
#'
#' @inheritParams pe_mediation
#'
#' @return A tidy `data.frame` of class `poem_tbl` with the same schema
#'   as [pe_mediation_linear()].
#'
#' @references
#' Guo, X., Li, R., Liu, J., and Zeng, M. (2024). Estimations and tests
#' for generalized mediation models with high-dimensional potential
#' mediators. \emph{Journal of Business & Economic Statistics, 42}(1),
#' 243--256. \doi{10.1080/07350015.2023.2174548}
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()], [pe_mediation_linear()],
#'   [pe_mediation_logistic()].
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
#'                              outcome = "count", c1 = 0.5, c2 = 0.4)
#' pe_mediation_poisson(d$X, d$Y, d$M)
#'
#' @export
pe_mediation_poisson <- function(X, Y, M, Z = NULL,
                                 method = c("Bonferroni", "BH", "BY"),
                                 scale = TRUE, error_level = 0.05,
                                 conf_level = 0.95, report_all_methods = FALSE,
                                 drop_constant = FALSE,
                                 lambda_grid = NULL,
                                 lambda_grid_reduced = NULL) {
  method <- match.arg(method)
  .validate_error_level(error_level)
  .validate_conf_level(conf_level)
  # report_all_methods = TRUE screens with all three multiplicity methods
  # (the chosen one stays primary) so their active sets can be compared.
  methods <- if (isTRUE(report_all_methods))
    unique(c(method, "Bonferroni", "BH", "BY")) else method
  inp <- .pe_validate_inputs(X, Y, M, Z, drop_constant = drop_constant)
  X <- inp$X; Y <- inp$Y; M <- inp$M; Z <- inp$Z
  n <- inp$n; q <- inp$q; s <- inp$s; p <- ncol(M)
  phi0 <- 1                          # known dispersion for the Poisson family
  lambda_grid <- .pe_resolve_grids(lambda_grid, lambda_grid_reduced, n, p,
                                   "count")$full

  # The outcome must be a nonnegative integer count with some variation.
  if (any(Y < 0, na.rm = TRUE))
    stop("`Y` must be nonnegative for the Poisson model.", call. = FALSE)
  if (any(abs(Y - round(Y)) > sqrt(.Machine$double.eps), na.rm = TRUE))
    stop("`Y` must contain integer counts for the Poisson model.", call. = FALSE)
  if (length(unique(Y[!is.na(Y)])) == 1L)
    stop("`Y` has no variation; the Poisson model is not identifiable.",
         call. = FALSE)

  if (!is.logical(scale) || length(scale) != 1L || is.na(scale))
    stop("`scale` must be TRUE or FALSE.", call. = FALSE)
  if (scale) {
    sc <- .pe_scale_data(X, Y, M, Z, center_y = FALSE)
    X <- sc$X; Y <- sc$Y; M <- sc$M; Z <- sc$Z
  }
  S <- Z

  # --- Penalized Poisson fit + HBIC tuning. (The Poisson tie-break takes the
  # first minimizer, matching the reference implementation.) ---
  ngrid <- length(lambda_grid)
  hbic <- numeric(ngrid)
  alpha_rcd <- matrix(0, nrow = ngrid, ncol = p + q + s)
  for (j in seq_len(ngrid)) {
    fit <- .HBIC_poisson(X, Y, M, lambda_grid[j], S)
    hbic[j] <- fit$HBIC
    alpha_rcd[j, ] <- fit$alpha_hat
  }
  id <- which(hbic == min(hbic))[1]
  alpha_hat <- alpha_rcd[id, ]
  alpha0_hat <- alpha_hat[1:p]
  alpha1_hat <- alpha_hat[(p + 1):(p + q)]
  alpha2_hat <- if (s > 0) alpha_hat[(p + q + 1):(p + q + s)] else NULL
  A <- which(alpha0_hat != 0)
  tuning <- .pe_tuning_record(lambda_grid, lambda_grid[id])

  if (length(A) == 0L)
    return(.pe_relabel_mediators(
      .pe_empty_output(n, p, q, "count (Poisson)", method,
                       error_level, conf_level, methods = methods,
                       tuning = tuning,
                       exposure_names = inp$exposure_names),
      inp$kept, inp$mediator_names))

  M_A <- as.matrix(M[, A])

  # --- Benchmark Wald inference for the Poisson model. ---
  hdmm <- .Testing_poisson(X, Y, M_A, S, phi0, alpha0_hat = alpha0_hat[A],
                           alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat)
  stat_hdmm <- as.numeric(hdmm$Sn)
  pval_hdmm <- .pe_chisq_pvalue(stat_hdmm, df = q)

  # --- Power enhancement component (outcome glm is Poisson). ---
  # --- Confidence interval, power enhancement screen(s), and tidy output. ---
  .pe_relabel_mediators(
    .pe_finalize(X, Y, M_A, S, A, stat_hdmm, pval_hdmm, hdmm$beta_hat,
                 hdmm$var_beta, p, n, q, "poisson", methods, error_level,
                 conf_level, "count (Poisson)", tuning = tuning,
                       exposure_names = inp$exposure_names),
    inp$kept, inp$mediator_names)
}
