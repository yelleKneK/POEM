# Estimation internals for the continuous-outcome (linear) PE test.
# Ported from the reference PEmediation implementation, preserving the exact
# numerics, with the reasoning behind each step documented. Not exported.

# Derivative of the SCAD penalty, evaluated for the local linear approximation
# (LLA). The folded-concave SCAD penalty is what gives the estimator its oracle
# property (it selects the true active set with probability approaching one and
# estimates the nonzero coefficients as if that set were known in advance),
# which is precisely the property the PE component relies on. The LLA replaces
# the nonconvex SCAD penalty by a weighted L1 penalty whose weight at each
# coefficient is the SCAD derivative there, so one lasso fit with these weights
# approximates the SCAD solution. With the usual SCAD shape parameter a = 3.7
# (Fan and Li, 2001), the derivative is 1 for small coefficients (full lasso
# shrinkage), tapers linearly to 0 over [lambda, a*lambda], and is 0 beyond
# (large coefficients are left unpenalized, removing lasso's bias).
#' @keywords internal
#' @noRd
.deSCAD <- function(z, lamb, a = 3.7) {
  1 * (z <= lamb) + pmax((a * lamb - z), 0) / ((a - 1) * lamb) * (lamb < z)
}

# One LLA pass for the linear mediator-outcome model: an initial lasso fit
# (Step 1) provides coefficients at which the SCAD derivatives are evaluated,
# and a second weighted lasso fit (Step 2) with those derivatives as penalty
# factors yields the LLA-SCAD estimate. The mediators M are penalized; the
# exposures X and confounders S (the last q + s columns) are not, so their
# penalty factors are left at 0. n_imp leaves the last n_imp mediator
# coefficients unpenalized (important mediators known a priori); it is 0 here.
# The fit is through the origin (intercept = FALSE) because the data are
# centered and standardized upstream.
#' @keywords internal
#' @noRd
.LLA_h1 <- function(X, Y, M, lamb, n, p, q, n_imp = 0, S = NULL, MV = NULL) {
  s <- if (length(S) == 0) 0L else ncol(S)
  # The design matrix [M, X, S] is identical for every lambda and for both LLA
  # steps, so build it once here (or reuse one passed in from the lambda loop)
  # rather than cbind()-ing it on each of the ~80 glmnet calls per fit. This is
  # numerically identical to the per-call cbind; it only avoids reallocation.
  if (is.null(MV)) MV <- if (s == 0L) cbind(M, X) else cbind(M, X, S)
  # Step 1: lasso (alpha = 1) penalizing only the mediator block.
  w1 <- matrix(0, nrow = (p + q + s), ncol = 1)
  w1[1:(p - n_imp)] <- 1
  alpha_int <- stats::coef(
    glmnet::glmnet(MV, Y, family = "gaussian", alpha = 1,
                   lambda = lamb, penalty.factor = w1, intercept = FALSE))[-1]
  # Step 2: SCAD derivatives at the Step 1 coefficients become the weights.
  w2 <- matrix(0, nrow = (p + q + s), ncol = 1)
  for (j in 1:(p - n_imp)) w2[j] <- .deSCAD(alpha_int[j], lamb)
  alpha <- stats::coef(
    glmnet::glmnet(MV, Y, family = "gaussian", alpha = 1,
                   lambda = lamb, penalty.factor = w2, intercept = FALSE))
  alpha[-1]
}

# High-dimensional BIC (HBIC) score at one lambda, used to choose the tuning
# parameter. Ordinary BIC underpenalizes model size when p grows with n; the
# HBIC of Wang, Kim, and Li (2013) replaces BIC's log(n) multiplier with
# log(log n) * log(p + q + s), which gives consistent selection in the
# high-dimensional regime. The score is log(residual variance) plus that
# size penalty times the model degrees of freedom (number of nonzero mediator
# coefficients, plus the always-included q exposures and s confounders). The
# lambda minimizing HBIC is selected.
#' @keywords internal
#' @noRd
.HBIC_calc_linear <- function(lamb, xx, yy, mm, S = NULL, n_imp = 0, MV = NULL) {
  n <- nrow(xx); p <- ncol(mm); q <- ncol(xx)
  if (is.null(S)) {
    s <- 0
    result <- .LLA_h1(xx, yy, mm, lamb, n, p, q, n_imp = n_imp, MV = MV)
    alpha0 <- result[1:p]; alpha1 <- result[(p + 1):(p + q)]; alpha2 <- NULL
    tmp <- yy - mm %*% alpha0 - xx %*% alpha1
  } else {
    s <- ncol(S)
    result <- .LLA_h1(xx, yy, mm, lamb, n, p, q, n_imp = n_imp, S = S, MV = MV)
    alpha0 <- result[1:p]; alpha1 <- result[(p + 1):(p + q)]
    alpha2 <- result[(p + q + 1):(p + q + s)]
    tmp <- yy - mm %*% alpha0 - xx %*% alpha1 - S %*% alpha2
  }
  df <- length(which(alpha0 != 0)) + q + s
  sigma_hat <- t(tmp) %*% tmp / n
  bic <- log(sigma_hat) + df * log(log(n)) * log(p + q + s) / n
  list(BIC = bic, alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
       sigma1_hat = sigma_hat)
}

# Wald inference for the total indirect effect in the linear mediation model.
# Given the penalized estimates restricted to the selected mediators, this
# computes the benchmark statistic S_n = n * beta' Sigma_beta^{-1} beta, the
# total indirect effect estimate beta = gamma_x - alpha_x (total effect minus
# direct effect), and the asymptotic covariance Sigma_beta. The covariance has
# two pieces (Guo et al., 2022): a sigma2 * Sigma_xx^{-1} part from the
# direct-path noise and a sigma1 * B part propagated through the X -> M map,
# where B is the inflation from regressing the mediators out. The confounder
# branch (s > 0) carries the same algebra through with V = [X, S] in place of
# X. S_n converges to chi-square_q under the no-total-indirect-effect null.
#' @keywords internal
#' @noRd
.inference_linear <- function(X, Y, M, S = NULL, M_tld = NULL,
                              alpha0_hat = NULL, alpha0_tld = NULL,
                              alpha1_hat = NULL, alpha2_hat = NULL,
                              alpha2_tld = NULL) {
  p <- ncol(M); q <- ncol(X); n <- nrow(X)
  if (is.null(S)) {
    # No confounders.
    Z <- cbind(M, X); s <- 0
    MS <- if (is.null(M_tld)) M else M_tld
    if (is.null(alpha0_tld)) alpha0_tld <- solve(t(MS) %*% MS) %*% t(MS) %*% Y
    RSS02 <- t(Y - MS %*% alpha0_tld) %*% (Y - MS %*% alpha0_tld)
  } else {
    Z <- cbind(M, X, S); s <- ncol(S)
    MS <- if (is.null(M_tld)) cbind(M, S) else cbind(M_tld, S)
    if (is.null(alpha0_tld)) alpha0_tld <- solve(t(MS) %*% MS) %*% t(MS) %*% Y
    RSS02 <- t(Y - MS %*% c(alpha0_tld, alpha2_tld)) %*%
             (Y - MS %*% c(alpha0_tld, alpha2_tld))
  }
  if (is.null(alpha0_hat)) {
    alpha_rf <- solve(t(Z) %*% Z) %*% t(Z) %*% Y
    alpha0_hat <- alpha_rf[1:p]; alpha1_hat <- alpha_rf[(p + 1):(p + q)]
    if (s > 0) alpha2_hat <- alpha_rf[(p + q + 1):(p + q + s)]
  } else {
    alpha_rf <- c(alpha0_hat, alpha1_hat, alpha2_hat)
  }

  res <- Y - Z %*% alpha_rf
  RSS12 <- as.numeric(t(res) %*% res)         # residual SS of the full model
  df <- p + q + s
  sigma1_hat <- RSS12 / (n - df)
  Sigma_MM <- t(M) %*% M / n

  if (s == 0) {
    gamma_hat <- solve(t(X) %*% X) %*% t(X) %*% Y
    sigmaT_hat <- t(Y - X %*% gamma_hat) %*% (Y - X %*% gamma_hat) / (n - q)
    beta_hat <- gamma_hat - alpha1_hat        # total minus direct = indirect
    sigma2_hat <- pmax(0, (sigmaT_hat - sigma1_hat))
    invXX <- solve(t(X) %*% X / n)
    Sigma_MX <- t(M) %*% X / n
    B <- invXX %*% t(Sigma_MX) %*%
         solve(Sigma_MM - Sigma_MX %*% invXX %*% t(Sigma_MX)) %*%
         Sigma_MX %*% invXX
    var_alpha1_hat <- sigma1_hat * (invXX + B)
    cov_beta_hat <- sigma2_hat * invXX + sigma1_hat * B
  } else {
    V <- cbind(X, S)
    gamma_hat <- solve(t(V) %*% V) %*% t(V) %*% Y
    sigmaT_hat <- t(Y - V %*% gamma_hat) %*% (Y - V %*% gamma_hat) / (n - q - s)
    beta_hat <- gamma_hat[1:q] - alpha1_hat
    sigma2_hat <- pmax(0, (sigmaT_hat - sigma1_hat))
    invVV <- solve(t(V) %*% V / n)
    Sigma_MV <- t(M) %*% V / n
    Sigma_VM <- t(Sigma_MV)
    B <- invVV %*% Sigma_VM %*%
         solve(Sigma_MM - Sigma_MV %*% invVV %*% Sigma_VM) %*%
         Sigma_MV %*% invVV
    var_alpha1_hat <- sigma1_hat * (invVV + B)[1:q, 1:q]
    cov_beta_hat <- sigma2_hat * invVV + sigma1_hat * B
    cov_beta_hat <- cov_beta_hat[1:q, 1:q]
  }

  Sn <- n * t(beta_hat) %*% solve(cov_beta_hat) %*% beta_hat   # Wald statistic
  Tn1 <- (n - df) * (RSS02 - RSS12) / RSS12                    # direct-effect LRT
  list(Sn = Sn, Tn = Tn1, beta_hat = beta_hat, alpha0_hat = alpha0_hat,
       alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat, B = B,
       var_beta = cov_beta_hat, var_alpha1_hat = var_alpha1_hat)
}
