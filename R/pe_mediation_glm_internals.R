# Estimation internals for the binary (logistic) and count (Poisson) PE
# tests, i.e. the generalized-outcome models. These routines are ported
# verbatim from the reference PEmediation implementation so the numerics match
# exactly; only the function names are changed to dotted internal names. They
# are documented as a group here rather than re-commented line by line, since
# the math mirrors the linear case (a penalized partial-likelihood fit, an HBIC
# tuning score, and a Wald inference for the total indirect effect) with the
# identity link replaced by the logit or log link.
#
# Why a hand-written coordinate descent (.MedReg) for the logistic fit rather
# than glmnet: the partial penalized LIKELIHOOD estimator of Guo et al. (2024)
# leaves the exposures and confounders unpenalized while applying a folded
# concave (SCAD/MCP) penalty to the mediators, and needs the SCAD/MCP solution
# (not the lasso) so the oracle property the PE component relies on holds.
# .MedReg is an iteratively reweighted coordinate descent over the strong /
# eligible active sets (the standard glmnet-style organization) that delivers
# exactly that estimator. The Poisson fit instead uses a glmnet lasso warm
# start followed by an IRLS quadratic approximation solved by .ZouAlgo_h1 (the
# weighted-LLA step), which is the count-outcome analogue.
#
# Group members:
#   .SCAD, .MCP            coordinate-descent thresholding for the two penalties
#   .p_binomial, .dedeb    logistic mean and its variance function b (weights)
#   .MedReg                penalized IRLS coordinate descent (logistic fit)
#   .HBIC_bino             HBIC tuning score for the logistic fit
#   .Testing_bino          Wald inference and S_n for the logistic model
#   .Ln_binomial           logistic partial log-likelihood (per observation)
#   .dedeb_poi             Poisson variance function b = exp(eta)
#   .ZouAlgo_h1            weighted-LLA solver used inside the Poisson fit
#   .poisson_h1            one penalized Poisson fit at a given lambda
#   .HBIC_poisson          HBIC tuning score for the Poisson fit
#   .Testing_poisson       Wald inference and S_n for the Poisson model
#   .Ln_poisson            Poisson partial log-likelihood (per observation)
#
# (.deSCAD, the SCAD derivative used by the Poisson LLA step, is defined once in
# pe_mediation_linear_internals.R and reused here.)
#
# @keywords internal
# @noRd

# ============================ LOGISTIC ============================
.SCAD <- function(z, l1, l2, gamma,v) {
  if (abs(z) <= l1)
  {return(0)}
  else if (abs(z) <= (l1*(1+l2)+l1))
  {return(sign(z)*(abs(z)-l1)/(v*(1+l2)))}
  else if (abs(z) <= gamma*l1*(1+l2))
  {return(sign(z)*(abs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2)))}
  else
  {return(z/(v*(1+l2)))}
}

# Internal helper:
#' @keywords internal
#' @noRd
.MCP <- function(z, l1, l2, gamma, v){
if (abs(z) <= l1) return(0)
else if (abs(z) <= gamma*l1*(1+l2)) return(sign(z)*(abs(z)-l1)/(v*(1+l2-1/gamma)))
else return(z/(v*(1+l2)))
}


# Internal helper:
#' @keywords internal
#' @noRd
.p_binomial <- function(eta) {
  if (eta > 10) {
    return(1)
  } else if (eta < -10) {
    return(0)
  } else {
    return(exp(eta)/(1+exp(eta)))
  }
}

# Internal helper:
#' @keywords internal
#' @noRd
.dedeb <- function(z){
  return(diag(c(exp(z)/(1+exp(z))^2)))
}

# Internal helper:
#' @keywords internal
#' @noRd
.MedReg <- function(X, y, family=c("gaussian","binomial","poisson"), penalty=c("MCP", "SCAD", "lasso"),
                   gamma=switch(penalty, SCAD=3.7, 3), alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100,
                   lambda, eps=1e-4, max_iter=10000, convex=TRUE, dfmax=p+1, penalty.factor=rep(1, ncol(x)),
                   warn=TRUE, returnX, Intercept = FALSE){
  # .MedReg fits the penalized generalized linear model by cyclic coordinate
  # descent with an active-set ("strong rule") strategy, the same organization
  # glmnet uses. The flow is: standardize the predictors; initialize at the
  # null (intercept-only) model; then for each tuning value lambda, repeatedly
  # (a) form the IRLS quadratic approximation of the log-likelihood at the
  # current fit (working response r and weights w), (b) cycle over the
  # currently active coordinates updating each by a soft/firm-threshold step
  # (.SCAD / .MCP), and (c) check whether any currently-inactive coordinate
  # should enter the active set, repeating until nothing new enters. The
  # mediator columns are penalized; the exposure/confounder columns carry a
  # zero penalty factor and are always retained.

  # Standardize predictors (center, scale to unit norm); std() records the
  # centers/scales so the coefficients can be returned on the original scale.
  x = ncvreg::std(X)
  p = ncol(x)
  n = nrow(x)
  ns = attr(x, "nonsingular")
  # Coersion
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  #if (class(x) != "matrix") {
  #  tmp <- try(x <- model.matrix(~0+., data=x), silent=TRUE)
  #  if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  #}
  if (storage.mode(x)=="integer") storage.mode(x) <- "double"
  if (!is.numeric(y)) {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must numeric or able to be coerced to numeric")
  }
  if (storage.mode(penalty.factor) != "double") storage.mode(penalty.factor) <- "double"
  # Preallocate the coefficient paths and the IRLS working vectors: b holds the
  # coefficient for every (lambda, predictor); b0 the intercepts; Dev the
  # deviance; e1/e2 the active-set flags; eta the linear predictor; r, s, w the
  # IRLS working response, residual, and weights; z the gradient used by the
  # strong rule.
  L = length(lambda)
  a = matrix(0,nrow = 1,ncol = p)
  b = matrix(0,nrow = L,ncol = p)
  lam = lambda
  b0 = matrix(0,nrow = 1,ncol = L)
  Dev = matrix(0, nrow = 1, ncol = L)
  e1 = matrix(0,nrow = 1, ncol = p)
  e2 = e1
  Eta = matrix(0,nrow = n, ncol = L)
  eta = matrix(0,nrow = n,ncol = 1)
  iter = matrix(0,nrow = 1,ncol = L)
  tot_iter = 0
  r = matrix(0,nrow = n,ncol = 1)
  s = matrix(0,nrow = n,ncol = 1)
  w = matrix(0,nrow = n,ncol = 1)
  z = matrix(0,nrow = 1,ncol = p)
  user = as.integer(any(penalty.factor==0))
  if (family=="binomial") {
    ## Initialization: start at the intercept-only logistic model, where the
    ## fitted probability is the sample mean ybar and the intercept is its
    ## log-odds. nullDev is the deviance of that null model (the yardstick for
    ## the saturation check below); s is the score (raw residual); z is the
    ## gradient of the log-likelihood at the null, used by the strong rule to
    ## decide which coordinates are eligible to enter the model.
    ybar = mean(y)
    b0[1] = log(ybar/(1-ybar))
    a0 = b0[1]
    nullDev = sum( -y*log(ybar) - (1-y)*log(1-ybar))
    s = y - ybar
    eta = rep(a0,n)
    z = t(s) %*% x

    ## If lam[0]=lam_max, skip lam[0]; closed form sol'n available
    if (user) {
      lstart = 1
    } else {
      lstart = 2
      Dev[1] = nullDev
      Eta[1,] = eta
    }
    for (l in lstart:L) {
      if (l != 1) {
        ## Assign a, a0
        a0 = b0[l-1]
        a =  b[(l-1),]

        ## Check dfmax
        nv = length(which(a!=0))
        if ((nv > dfmax) | (tot_iter == max_iter)) {
          for (ll in l:L) iter[ll] = NA
          break
        }

        ## Determine eligible set
        if (penalty== "MCP") cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lam[l-1])
        if (penalty== "SCAD") cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lam[l-1])
        for (j in 1:p) if (abs(z[j]) > (cutoff * alpha * penalty.factor[j])) e2[j] = 1
      } else {

        ## Determine eligible set
        lmax = max(abs(z))
        if (penalty == "MCP") cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lmax)
        if (penalty == "SCAD") cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lmax)
        for (j in 1:p){
          if (abs(z[j]) > (cutoff * alpha *  penalty.factor[j])){
            e2[j] = 1
          }
        }
      }


      while (tot_iter < max_iter) {
        while (tot_iter < max_iter) {
          maxChange = 1e-3
          while ((tot_iter < max_iter) & (maxChange > eps)) {
            iter[l] = iter[l] + 1
            tot_iter = tot_iter + 1
            Dev[l] = 0
            # IRLS reweighting: at the current linear predictor eta, compute the
            # fitted probability pi, the binomial working weight w = pi(1-pi)
            # (floored away from 0 for numerical stability), the score s, and
            # the working response r = s/w. Accumulate the deviance Dev as we go.
            for (i in 1:n) {
              pi = .p_binomial(eta[i])
              w[i] = pmax(pi*(1-pi), 0.0001)
              s[i] = y[i] - pi
              r[i] = s[i]/w[i]
              if (y[i]==1) Dev[l] = Dev[l] - log(pi)
              if (y[i]==0) Dev[l] = Dev[l] - log(1-pi)
            }
            if (Dev[l]/nullDev < 0.01) {
              if (warn) warning("Model saturated; exiting...");
              tot_iter = max_iter
              break
            }

            ## Intercept
            if(Intercept){
              xwr = t(w)%*%r
              xwx = sum(w)
              b0[l] = xwr/xwx + a0
              for (i in 1:n) {
                si = b0[l] - a0
                r[i] =  r[i] -si
                eta[i] = eta[i]+ si
              }
              maxChange = abs(si)*xwx/n
            }else{
              maxChange = 0
            }
            ## Covariates: cycle over the active coordinates (e1[j] == 1). For
            ## each, form the weighted single-coordinate least-squares target u
            ## (the gradient plus the current contribution) and curvature v,
            ## then take the firm-thresholding step for the chosen penalty
            ## (.MCP or .SCAD); l1/l2 are the lasso/ridge parts of the elastic
            ## net penalty. The coefficient that is unpenalized (penalty.factor
            ## 0, the exposures and confounders) gets l1 = 0 and so is updated
            ## by ordinary least squares.
            for (j in 1:p) {
              if (e1[j]) {

                ## Calculate u, v
                xwr = sum(x[,j] * r * w)
                xwx = sum(x[,j]^2*w)
                u = xwr/n + (xwx/n)*a[j]
                v = xwx/n

                ## Update b_j
                l1 = lam[l] * penalty.factor[j] * alpha
                l2 = lam[l] * penalty.factor[j] * (1-alpha)
                if (penalty == "MCP") b[l,j] = .MCP(u, l1, l2, gamma, v)
                if (penalty =="SCAD") b[l,j] = .SCAD(u, l1, l2, gamma, v)

                ## Update r: propagate the coefficient change through the
                ## residual and linear predictor, and track the largest change
                ## for the convergence test.
                shift = b[l,j] - a[j]
                if (shift !=0) {
                  for (i in 1:n) {
                    si = shift*x[i,j]
                    r[i] =r[i]- si
                    eta[i] =eta[i]+ si
                  }
                  if (abs(shift)*sqrt(v) > maxChange) {
                    maxChange = abs(shift)*sqrt(v)}
                }
              }
            }

            ## Check for convergence
            a0 = b0[l]
            a = b[l,]

          }
          ## Scan for violations in strong set
          violations = 0
          for (j in 1:p) {
            if (e1[j]==0 & e2[j]==1) {
              z[j] = t(x[,j])%*%s/n
              l1 = lam[l] * penalty.factor[j] * alpha
              if (abs(z[j]) > l1) {
                e1[j] = 1
                e2[j] = 1
                violations = violations + 1
              }
            }
          }
          if (violations==0) break
        }

        ## Scan for violations in rest
        violations = 0
        for (j in 1:p) {
          if (e2[j]==0) {
            z[j] = t(x[,j])%*%s/n
            l1 = lam[l] * penalty.factor[j] * alpha
            if (abs(z[j]) > l1) {
              e1[j] = 1
              e2[j] = 1
              violations = violations+1
            }
          }
        }
        if (violations==0) {
          Eta[,l] = eta
          break
        }
      }
    }
    ## Unstandardize
    b =matrix(b, p, L,byrow=TRUE)

    if(Intercept){
      beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda))
      bb <- b/attr(x, "scale")[ns]
      beta[ns+1,] <- bb
      beta[1,] <- b0 - attr(x, "center")[ns]%*% bb
    }else{
      beta = b/attr(x, "scale")
    }
  }
  return(list(b = b,beta = beta,Dev=Dev, Eta = Eta, iter = iter))
}


# Internal helper:
#' @keywords internal
#' @noRd
.HBIC_bino <- function(X,Y,M,S = NULL,w,lamb,
                      Y_family = "binomial",
                      penalty_type = "SCAD"){

  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(M)

  if (is.null(S)) {
    s <- 0
    x <- as.matrix(cbind(M, X))
  } else {
    S <- as.matrix(S)
    if (nrow(S) != n) stop("S must have the same number of rows as X.")
    s <- ncol(S)
    x <- as.matrix(cbind(M, X, S))
  }

  #x = as.matrix(cbind(M,X))
  #p = ncol(x)
  #n = nrow(x)

  result <- .MedReg(
    x, Y,
    family = Y_family,
    penalty = penalty_type,
    lambda = lamb,
    penalty.factor = w,
    alpha = 1,
    Intercept = FALSE
  )

  alpha <- as.vector(result$beta)

  alpha0 <- alpha[1:p]
  alpha1 <- alpha[(p + 1):(p + q)]
  alpha2 <- if (s > 0) alpha[(p + q + 1):(p + q + s)] else NULL

  logLkhd <- .Ln_binomial(X, Y, M, alpha0, alpha1, S = S, alpha2 = alpha2)

  # mediator coefficients are sparse; X and S are always counted
  df <- sum(abs(alpha0) > 0) + q + s

  # p + q + s is the total number of regression coefficients
  bic <- -n * logLkhd + df * log(log(n)) * log(p + q + s)

  return(list(
    HBIC = bic,
    logLkhd = logLkhd,
    alpha = alpha,
    alpha0 = alpha0,
    alpha1 = alpha1,
    alpha2 = alpha2
  ))
}


# Internal helper:
#' @keywords internal
#' @noRd
.Testing_bino <-function(X,Y,M_A, S = NULL, phi0=1,
                        a0_hat=NULL, a1_hat=NULL,a2_hat = NULL,
                        Y_family = "binomial"){

  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(M_A)

  if (!is.null(S)) {
    s <- ncol(S)
  } else {
    s <- 0
  }

  # if(is.null(a0_hat)){
  # #Ocl = glm(Y~0 + M_A + X,family = binomial(link = "logit")) # Refit the model with selected subset of variable
  # Ocl = glm(Y~0 + M_A + X,family = Y_family) # Refit the model with selected subset of variable
  # a0 = Ocl$coefficients[1:p]
  # a1 = Ocl$coefficients[(p+1):(p+q)]
  # }
  # else{a0 = a0_hat; a1 = a1_hat}

  # Refit logistic model if coefficients are not supplied
  if (is.null(a0_hat) || is.null(a1_hat) || (!is.null(S) && is.null(a2_hat))) {
    dat <- data.frame(Y = Y, M_A, X)
    if (!is.null(S)) dat <- cbind(dat, S)

    # no intercept, consistent with your original setup
    Ocl <- glm(Y ~ 0 + ., family = Y_family, data = dat)

    cf <- coef(Ocl)
    a0 <- cf[1:p]
    a1 <- cf[(p + 1):(p + q)]
    a2 <- if (!is.null(S)) cf[(p + q + 1):(p + q + s)] else NULL
  } else {
    a0 <- a0_hat
    a1 <- a1_hat
    a2 <- if (!is.null(S)) a2_hat else NULL
  }


  if(is.null(S))
  {
    # Total indirect effect on the logit scale, beta = Gamma_x alpha_m,
    # estimated by regressing the fitted mediator signal M_A %*% a0 on the
    # exposure; hatsigma2 is the residual variance of that regression (the
    # X -> M noise that feeds the delta-method variance below).
    beta_hat = solve(t(X)%*%X) %*% t(X) %*% M_A %*% a0
    hatsigma2 = sum((M_A %*% a0 - X %*% beta_hat)^2)/(n-1)
    # Fisher-information blocks of the logistic fit: bbb = b''(eta) = pi(1-pi)
    # is the GLM variance function, and the Sigma* matrices are the
    # information-weighted cross products. Sigmamm.x is the mediator block after
    # partialling out the exposure; the two-term var_beta combines the X -> M
    # noise (hatsigma2) with the M -> Y estimation uncertainty (phi0 term).
    z = X %*% a1 + M_A %*% a0
    bbb = .dedeb(z) # the second derivative of b(z)
    Sigmaxx = t(X) %*% bbb %*% X/n
    Sigmaxm = t(X) %*% bbb %*% M_A/n
    Sigmamm = t(M_A) %*% bbb %*% M_A/n
    tldSigmaxx = t(X)%*%X/n
    tldSigmaxm = t(X) %*%M_A/n
    Sigmaxx_inv = solve(Sigmaxx)
    Sigmamm.x = Sigmamm - t(Sigmaxm) %*% Sigmaxx_inv %*% Sigmaxm
    Sigmamm.x_inv = solve(Sigmamm.x)
    tldSigmaxx_inv = solve(tldSigmaxx)
    # Calculate the variance of \hat{alpha_1} and \hat{\beta_1}
    var_alpha = phi0*(Sigmaxx_inv +Sigmaxx_inv %*%Sigmaxm %*% Sigmamm.x_inv %*% t(Sigmaxm) %*%Sigmaxx_inv)
    var_beta = hatsigma2 * tldSigmaxx_inv + phi0* tldSigmaxx_inv %*% tldSigmaxm %*% Sigmamm.x_inv %*% t(tldSigmaxm) %*%tldSigmaxx_inv

  } else {

    # Residual-maker for S
    P_S <- S %*% solve(t(S) %*% S) %*% t(S)
    M_S <- diag(n) - P_S

    # beta_hat: coefficient of X in regression of M_A %*% a0 on (X, S)
    beta_hat <- solve(t(X) %*% M_S %*% X) %*% t(X) %*% M_S %*% M_A %*% a0

    # coefficient of S in same regression
    gamma_hat <- solve(t(S) %*% S) %*% t(S) %*% (M_A %*% a0 - X %*% beta_hat)

    # residual variance
    resid_vec <- M_A %*% a0 - X %*% beta_hat - S %*% gamma_hat
    hatsigma2 <- sum(resid_vec^2) / (n - ncol(X) - ncol(S))

    # logistic linear predictor
    z <- X %*% a1 + M_A %*% a0 + S %*% a2

    # second derivative matrix
    bbb <- .dedeb(z)

    # weighted blocks
    Sigmaxx <- t(X)   %*% bbb %*% X   / n
    Sigmaxs <- t(X)   %*% bbb %*% S   / n
    Sigmaxm <- t(X)   %*% bbb %*% M_A / n

    Sigmass <- t(S)   %*% bbb %*% S   / n
    Sigmasm <- t(S)   %*% bbb %*% M_A / n

    Sigmamm <- t(M_A) %*% bbb %*% M_A / n

    # unweighted partialed-out matrices
    tldSigmaxx.s <- t(X) %*% M_S %*% X   / n
    tldSigmaxm.s <- t(X) %*% M_S %*% M_A / n
    tldSigmaxx.s_inv <- solve(tldSigmaxx.s)

    # block matrices for (X, S)
    Sigma_xs_xs <- rbind(
      cbind(Sigmaxx, Sigmaxs),
      cbind(t(Sigmaxs), Sigmass)
    )

    Sigma_xs_m <- rbind(
      Sigmaxm,
      Sigmasm
    )

    # mediator block adjusted for X and S
    Sigmamm.xs <- Sigmamm - t(Sigma_xs_m) %*% solve(Sigma_xs_xs) %*% Sigma_xs_m
    Sigmamm.xs_inv <- solve(Sigmamm.xs)

    # X block adjusted for S
    Sigmaxx.s <- Sigmaxx - Sigmaxs %*% solve(Sigmass) %*% t(Sigmaxs)
    Sigmaxm.s <- Sigmaxm - Sigmaxs %*% solve(Sigmass) %*% Sigmasm
    Sigmaxx.s_inv <- solve(Sigmaxx.s)

    # variance of alpha1
    var_alpha <- phi0 * (
      Sigmaxx.s_inv +
        Sigmaxx.s_inv %*% Sigmaxm.s %*% Sigmamm.xs_inv %*%
        t(Sigmaxm.s) %*% Sigmaxx.s_inv
    )

    # variance of beta_hat
    var_beta <- hatsigma2 * tldSigmaxx.s_inv +
      phi0 * tldSigmaxx.s_inv %*% tldSigmaxm.s %*% Sigmamm.xs_inv %*%
      t(tldSigmaxm.s) %*% tldSigmaxx.s_inv
  }


  Sn <- n * t(beta_hat) %*% solve(var_beta) %*% beta_hat

  # constrained model: drop X, keep M_A and S
  if (is.null(S)) {
    refit_tld <- glm(Y ~ 0 + M_A, family = Y_family)
    alpha0_tld <- as.matrix(coef(refit_tld))
    alpha2_tld <- NULL
    Tn <- 2 * n * (.Ln_binomial(X, Y, M_A, a0, a1) -
                     .Ln_binomial(X, Y, M_A, alpha0_tld, as.matrix(integer(q)))) / phi0
  } else {
    dat_tld <- data.frame(Y = Y, M_A, S)
    refit_tld <- glm(Y ~ 0 + ., family = Y_family, data = dat_tld)
    cf_tld <- coef(refit_tld)
    alpha0_tld <- as.matrix(cf_tld[1:p])
    alpha2_tld <- as.matrix(cf_tld[(p + 1):(p + s)])

    Tn <- 2 * n * (
      .Ln_binomial(X, Y, M_A, a0, a1, S = S, alpha2 = a2) -
        .Ln_binomial(X, Y, M_A, alpha0_tld, as.matrix(integer(q)), S = S, alpha2 = alpha2_tld)
    ) / phi0
  }

  # # Test of direct effect likelihood ratio test
  # Tn = 2*n*(.Ln_binomial(X,Y,M_A,a0,a1) - .Ln_binomial(X,Y,M_A,alpha0_tld,as.matrix(integer(q))))/phi0
  # #Tn2 = 2 * n * (logLik(Ocl) - logLik(refit_tld))/phi0
  # return(list(beta_hat = beta_hat, alpha0_hat = a0, alpha1_hat = a1,
  #             var_beta = var_beta, var_alpha1 = var_alpha,
  #             Sn = Sn, Tn = Tn, sigma2 = hatsigma2))
  #

  return(list(
    beta_hat = beta_hat,
    alpha0_hat = a0,
    alpha1_hat = a1,
    alpha2_hat = a2,
    var_beta = var_beta,
    var_alpha1 = var_alpha,
    Sn = Sn,
    Tn = Tn,
    sigma2 = hatsigma2
  ))
}

# Internal helper:
#' @keywords internal
#' @noRd
.Ln_binomial <- function(X, Y, M, alpha0, alpha1, S = NULL, alpha2 = NULL) {
  eta <- (M%*%alpha0 + X%*% alpha1)
  if (!is.null(S)) {
    eta <- eta + S %*% alpha2
  }
  return(mean(diag(Y)%*% eta - log(1 + exp(eta))))
}

# ============================= POISSON ============================
.dedeb_poi<-function(z){
  return(exp(z))
}


# Internal helper:
#' @keywords internal
#' @noRd
.ZouAlgo_h1 <- function(X,Y,M,w,lamb){
  X = as.matrix(X) # unpenalized block: can be X or cbind(X, S)
  M = as.matrix(M) # penalized mediator block
  n = nrow(X)
  p = ncol(M)
  q = ncol(X)
  alpha_int = matrix(NA,ncol=1,nrow=(p+q))
  U = which(w == 0)
  V = which(w!=0)
  Xt = sqrt(2)*cbind(M,X)
  Xts = Xt
  for(j in 1:length(V)){
    Xts[,V[j]] = Xt[,V[j]] * lamb/w[V[j]]
  }
  Xus = as.matrix(Xts[,U])
  Xvs = as.matrix(Xts[,V])
  Pu = Xus%*%solve(t(Xus)%*%Xus)%*%t(Xus)
  Qu = diag(n) - Pu
  Ys = sqrt(2)*Y
  Yss1 = sqrt(2)*Qu%*%Y
  Xvss1 = Qu%*%Xvs
  reg1 = glmnet(Xvss1,Yss1,family = "gaussian",alpha=1,lambda=lamb)
  Betavs = matrix(reg1$beta,ncol=1)
  Betaus = solve(t(Xus)%*%Xus)%*%t(Xus)%*%(Ys - Xvs%*%Betavs)
  alpha_int[U] = Betaus
  alpha_int[V] = Betavs*lamb/w[V]
  return(alpha_int)
}



# Internal helper:
#' @keywords internal
#' @noRd
.poisson_h1 <- function(X,Y,M,lamb,S=NULL){
  p = ncol(M)
  q = ncol(X)
  n = nrow(M)

  if(is.null(S)){
    s = 0
    Z = as.matrix(cbind(M,X))
  }else{
    s = ncol(S)
    Z = as.matrix(cbind(M,X,S))
  }

  # Step 1 using Lasso
  # penalize M, do not penalize X or S
  w1 = vector(mode = "double",length= (p+q+s))
  w1[1:p] = 1
  res1 = glmnet(Z, Y,family = "poisson", alpha=1, lambda=lamb,penalty.factor = w1)
  alpha_int = as.matrix(res1$beta)

  # Step 2 using linear approximation of SCAD
  w2 =vector(mode = "double",length= (p+q+s))
  for(j in 1:p){
    w2[j] = .deSCAD(alpha_int[j],lamb)
  }
  xs = Z
  ys = Y
  for (i in 1:nrow(M)){
    b_dede = sqrt(as.numeric(exp(t(Z[i,])%*%alpha_int)))
    xs[i,] = b_dede%*%Z[i,]
    ys[i] = b_dede* t(Z[i,])%*% alpha_int + (Y[i] - exp(t(Z[i,])%*%alpha_int))/b_dede
  }
  alpha = .ZouAlgo_h1(xs[,(p+1):(p+q+s)],ys, xs[,1:p],w2,lamb)
  return(alpha)
}

# Internal helper:
#' @keywords internal
#' @noRd
.HBIC_poisson <-function(X,Y,M,lamb,S = NULL){
  M = as.matrix(M)
  X = as.matrix(X)
  p = ncol(M)
  q = ncol(X)
  n = nrow(M)
  #Cn = log(n*log(n))/2
  if(is.null(S))
  {
    s <- 0
    Z = as.matrix(cbind(M,X))
    alpha_hat = .poisson_h1(X,Y,M,lamb)
  }else{
    s <- ncol(S)
    Z <- as.matrix(cbind(M, X, S))
    alpha_hat <- .poisson_h1(X, Y, M, lamb, S = S)
  }

  eta <- Z%*%alpha_hat
  ETA_MAX <- log(.Machine$double.xmax) - 1
  eta_safe <- pmin(eta, ETA_MAX)
  logLkhd <- mean(Y*eta_safe - exp(eta_safe))

  tol <- 10e-12
  alpha0_hat <- alpha_hat[1:p, , drop = FALSE]
  df <- sum(abs(alpha0_hat) > tol) + q + s

  Cn = log(log(n))

  HBIC <- -logLkhd + df * Cn * log(p + q + s) / n

  return(list(HBIC=HBIC,alpha_hat=alpha_hat))
}

# Internal helper:
#' @keywords internal
#' @noRd
.Testing_poisson <- function(X,Y,M_A,S=NULL, phi0=1,
                            alpha0_hat=NULL, alpha1_hat=NULL, alpha2_hat=NULL){
  Y_family = "poisson"
  M_A = as.matrix(M_A)
  pm = ncol(M_A)
  q = ncol(X)
  n = nrow(M_A)

  if (!is.null(S)) {
    s <- ncol(S)
  } else {
    s <- 0
  }

  #if(is.null(alpha0_hat)){
  #refit = glm(Y ~ 0 + M_A+X,family = Y_family)
  #alpha0_hat = as.matrix(coef(refit)[1:pm])
  #alpha1_hat = as.matrix(coef(refit)[(pm+1):(pm+q)])
  #}
  #else{alpha0_hat = alpha0_hat; alpha1_hat = alpha1_hat}

  ## Refit if coefficient estimates are not supplied
  if (is.null(alpha0_hat) || is.null(alpha1_hat) || (!is.null(S) && is.null(alpha2_hat))) {
    if (is.null(S)) {
      refit <- glm(Y ~ 0 + M_A + X, family = Y_family)
      cf <- coef(refit)
      alpha0_hat <- as.matrix(cf[1:pm])
      alpha1_hat <- as.matrix(cf[(pm + 1):(pm + q)])
      alpha2_hat <- NULL
    } else {
      refit <- glm(Y ~ 0 + M_A + X + S, family = Y_family)
      cf <- coef(refit)
      alpha0_hat <- as.matrix(cf[1:pm])
      alpha1_hat <- as.matrix(cf[(pm + 1):(pm + q)])
      alpha2_hat <- as.matrix(cf[(pm + q + 1):(pm + q + s)])
    }
  } else {
    alpha0_hat <- as.matrix(alpha0_hat)
    alpha1_hat <- as.matrix(alpha1_hat)
    if (!is.null(S)) alpha2_hat <- as.matrix(alpha2_hat)
  }


  if(is.null(S))
  {

    beta_hat = solve(t(X)%*%X) %*% t(X) %*% M_A %*% alpha0_hat
    hatsigma2 = sum((M_A %*% alpha0_hat - X %*% beta_hat)^2)/(n-q)
    z = X %*% alpha1_hat + M_A %*% alpha0_hat
    bbb = diag(c(.dedeb_poi(z))) # the second derivative of b(z)
    Sigmaxx = t(X) %*% bbb %*% X/n
    Sigmaxm = t(X) %*% bbb %*% M_A/n
    Sigmamm = t(M_A) %*% bbb %*% M_A/n
    tldSigmaxx = t(X)%*%X/n
    tldSigmaxm = t(X) %*%M_A/n
    Sigmaxx_inv = solve(Sigmaxx)
    Sigmamm.x = Sigmamm - t(Sigmaxm) %*% Sigmaxx_inv %*% Sigmaxm
    Sigmamm.x_inv = solve(Sigmamm.x)
    tldSigmaxx_inv = solve(tldSigmaxx)
    var_alpha = phi0*(Sigmaxx_inv +Sigmaxx_inv %*%Sigmaxm %*% Sigmamm.x_inv %*% t(Sigmaxm) %*%Sigmaxx_inv)
    var_beta = hatsigma2 * tldSigmaxx_inv + phi0* tldSigmaxx_inv %*% tldSigmaxm %*% Sigmamm.x_inv %*% t(tldSigmaxm) %*%tldSigmaxx_inv

  } else {
    P_S <- S %*% solve(t(S) %*% S) %*% t(S)
    M_S <- diag(n) - P_S

    beta_hat <- solve(t(X) %*% M_S %*% X) %*% t(X) %*% M_S %*% M_A %*% alpha0_hat

    gamma_hat <- solve(t(S) %*% S) %*% t(S) %*% (M_A %*% alpha0_hat - X %*% beta_hat)
    resid_vec <- M_A %*% alpha0_hat - X %*% beta_hat - S %*% gamma_hat
    hatsigma2 <- sum(resid_vec^2) / (n - q - s)

    z = X %*% alpha1_hat + M_A %*% alpha0_hat + S %*% alpha2_hat
    bbb = diag(c(.dedeb_poi(z)))

    tldSigmaxx <- t(X) %*% M_S %*% X / n
    tldSigmaxm <- t(X) %*% M_S %*% M_A / n
    z <- X %*% alpha1_hat + M_A %*% alpha0_hat + S %*% alpha2_hat
    Sigmaxx <- t(X) %*% bbb %*% X / n
    Sigmaxm <- t(X) %*% bbb %*% M_A / n
    Sigmamm <- t(M_A) %*% bbb %*% M_A / n
    tldSigmaxx_inv <- solve(tldSigmaxx)
    Sigmaxs <- t(X) %*% bbb %*% S / n
    Sigmass <- t(S) %*% bbb %*% S / n
    Sigmasm <- t(S) %*% bbb %*% M_A / n
    Sigma_xs_xs <- rbind(cbind(Sigmaxx, Sigmaxs),cbind(t(Sigmaxs), Sigmass))
    Sigma_xs_m <- rbind(Sigmaxm,Sigmasm)

    Sigmamm.xs <- Sigmamm - t(Sigma_xs_m) %*% solve(Sigma_xs_xs) %*% Sigma_xs_m
    Sigmamm.xs_inv <- solve(Sigmamm.xs)

    ## X block adjusted for S
    Sigmaxx.s <- Sigmaxx - Sigmaxs %*% solve(Sigmass) %*% t(Sigmaxs)
    Sigmaxm.s <- Sigmaxm - Sigmaxs %*% solve(Sigmass) %*% Sigmasm
    Sigmaxx.s_inv <- solve(Sigmaxx.s)

    var_alpha <- phi0 * (
      Sigmaxx.s_inv +
        Sigmaxx.s_inv %*% Sigmaxm.s %*% Sigmamm.xs_inv %*%
        t(Sigmaxm.s) %*% Sigmaxx.s_inv
    )

    var_beta <- hatsigma2 * tldSigmaxx_inv +
      phi0 * tldSigmaxx_inv %*% tldSigmaxm %*% Sigmamm.xs_inv %*%
      t(tldSigmaxm) %*% tldSigmaxx_inv
  }

  # Test on indirect effect
  Sn =n * t(beta_hat) %*% solve(var_beta) %*% beta_hat

  # Test on direct effect
  ## constrained model removes X but keeps M_A and S
  if (is.null(S)) {
    refit_tld = glm(Y ~ 0 + M_A,family = Y_family)
    alpha0_tld = as.matrix(coef(refit_tld))
    Tn = 2*n*(.Ln_poisson(X,Y,M_A,alpha0_hat,alpha1_hat) - .Ln_poisson(X,Y,M_A,alpha0_tld,as.matrix(integer(q))))/phi0
  } else {
    refit_tld <- glm(Y ~ 0 + M_A + S, family = Y_family)
    cf_tld <- coef(refit_tld)
    alpha0_tld <- as.matrix(cf_tld[1:pm])
    alpha2_tld <- as.matrix(cf_tld[(pm + 1):(pm + s)])
    Tn <- 2*n*(
      .Ln_poisson(X, Y, M_A, alpha0_hat, alpha1_hat, S = S, a2 = alpha2_hat) -
        .Ln_poisson(X, Y, M_A, alpha0_tld, as.matrix(integer(q)), S = S, a2 = alpha2_tld)
    ) / phi0
  }
  return(list(beta_hat = beta_hat, alpha0_hat = alpha0_hat, alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat,
              var_beta = var_beta, var_alpha1 = var_alpha, Sn = Sn, Tn = Tn, sigma2 = hatsigma2))
}

# Internal helper:
#' @keywords internal
#' @noRd
.Ln_poisson <- function(X,Y,M_A,a0,a1,S = NULL, a2= NULL){
  z = X %*% a1 + M_A %*% a0
  if (!is.null(S)) {
    z <- z + S %*% a2
  }
  return(mean(diag(c(Y))%*%z - exp(z)))
}
