#' Power-enhanced mediation test for a binary outcome
#'
#' Tests the global null hypothesis of no active mediator for a binary
#' (0/1) outcome modeled with logistic regression, reporting both the
#' benchmark Wald test on the total indirect effect and the
#' power-enhanced (PE) test. This is the binary-outcome worker behind
#' [pe_mediation()].
#'
#' @details
#' The outcome follows a logistic mediation model,
#' \eqn{\mathrm{logit}\,P(Y = 1 \mid M, X, Z) = \alpha_m' M + \alpha_x' X
#' + \alpha_z' Z}, with a linear mediator model as in
#' [pe_mediation_linear()]. The mediator coefficients are estimated by
#' partial penalized likelihood with a SCAD penalty (a coordinate-descent
#' fit that leaves the exposures and confounders unpenalized), the tuning
#' parameter chosen by the high-dimensional BIC. On the logit (link)
#' scale the total indirect effect \eqn{\beta = \Gamma_x \alpha_m} and
#' the power enhancement construction are exactly as in the continuous
#' case (Guo et al., 2024). The total indirect effect here is on the
#' log-odds scale, not the probability scale.
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
#'   [pe_mediation_poisson()].
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' # Contrasting mediation with a binary outcome: the benchmark Wald test
#' # does not reject, the PE test does and identifies an active mediator.
#' d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
#'                              outcome = "binary", c1 = 1, c2 = 1)
#' pe_mediation_logistic(d$X, d$Y, d$M)
#'
#' @export
pe_mediation_logistic <- function(X, Y, M, Z = NULL,
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
  phi0 <- 1                          # known dispersion for the binomial family
  lambda_grid <- .pe_resolve_grids(lambda_grid, lambda_grid_reduced, n, p,
                                   "binary")$full

  # The outcome must be binary and coded 0/1 for the logistic model.
  y_unique <- unique(Y[!is.na(Y)])
  if (length(y_unique) != 2L)
    stop("`Y` must be binary (two distinct values) for the logistic model.",
         call. = FALSE)
  if (!all(y_unique %in% c(0, 1)))
    stop("`Y` must be coded as 0/1 for the logistic model. Detected values: ",
         paste(y_unique, collapse = ", "), ".", call. = FALSE)

  if (!is.logical(scale) || length(scale) != 1L || is.na(scale))
    stop("`scale` must be TRUE or FALSE.", call. = FALSE)
  if (scale) {
    sc <- .pe_scale_data(X, Y, M, Z, center_y = FALSE)
    X <- sc$X; Y <- sc$Y; M <- sc$M; Z <- sc$Z
  }
  S <- Z

  # --- Penalized partial-likelihood fit + HBIC tuning. The penalty factor w
  # penalizes only the mediator block (first p columns); exposures and
  # confounders are always retained. ---
  w <- rep(0, p + q + s); w[1:p] <- 1
  ngrid <- length(lambda_grid)
  hbic <- numeric(ngrid)
  alpha_rcd <- matrix(0, nrow = ngrid, ncol = p + q + s)
  for (j in seq_len(ngrid)) {
    fit <- .HBIC_bino(X = X, Y = Y, M = M, S = S, w = w, lamb = lambda_grid[j])
    hbic[j] <- fit$HBIC
    alpha_rcd[j, ] <- fit$alpha
  }
  id <- utils::tail(which(hbic == min(hbic)), 1)
  alpha_hat <- alpha_rcd[id, ]
  alpha0_hat <- alpha_hat[1:p]
  alpha1_hat <- alpha_hat[(p + 1):(p + q)]
  alpha2_hat <- if (s > 0) alpha_hat[(p + q + 1):(p + q + s)] else NULL
  A <- which(alpha0_hat != 0)
  tuning <- .pe_tuning_record(lambda_grid, lambda_grid[id])

  if (length(A) == 0L)
    return(.pe_relabel_mediators(
      .pe_empty_output(n, p, q, "binary (logistic)", method,
                       error_level, conf_level, methods = methods,
                       tuning = tuning,
                       exposure_names = inp$exposure_names),
      inp$kept, inp$mediator_names))

  M_A <- as.matrix(M[, A])

  # --- Benchmark Wald inference for the logistic model. ---
  hdmm <- .Testing_bino(X, Y, M_A, S, phi0, a0_hat = alpha0_hat[A],
                        a1_hat = alpha1_hat, a2_hat = alpha2_hat)
  stat_hdmm <- as.numeric(hdmm$Sn)
  pval_hdmm <- .pe_chisq_pvalue(stat_hdmm, df = q)

  # --- Power enhancement component (outcome glm is binomial). ---
  # --- Confidence interval, power enhancement screen(s), and tidy output. ---
  .pe_relabel_mediators(
    .pe_finalize(X, Y, M_A, S, A, stat_hdmm, pval_hdmm, hdmm$beta_hat,
                 hdmm$var_beta, p, n, q, "binomial", methods, error_level,
                 conf_level, "binary (logistic)", tuning = tuning,
                       exposure_names = inp$exposure_names),
    inp$kept, inp$mediator_names)
}
