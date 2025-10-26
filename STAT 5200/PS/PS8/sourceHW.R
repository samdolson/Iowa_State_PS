metropforgamma <- function(dat, start, priorpars, jumpvars, B, M){
  # Metropolis for a one-sample Gamma(shape = alpha, rate = beta) model
  # with product prior: alpha ~ Uniform(0, A), beta ~ Gamma(gamma0, lambda0)
  #
  # dat        : vector of observed positive data (y_i > 0)
  # start      : c(alpha0, beta0)  -- starting values for (alpha, beta)
  # priorpars  : c(gamma0, lambda0, A)
  #              - gamma0, lambda0 are shape/rate of prior on beta
  #              - A is the upper bound for alpha's Uniform(0, A) prior
  # jumpvars   : c(valpha, vbeta)  -- proposal variances for random-walk jumps
  # B          : burn-in iterations
  # M          : number of kept Monte Carlo draws
  #
  # Notes on parameterization/statistics:
  # - Likelihood: Y_i ~ Gamma(alpha, beta) with density
  #       f(y | alpha, beta) = beta^alpha / Gamma(alpha) * y^(alpha-1) * exp(-beta*y)
  #   The log-likelihood is computed in a numerically stable way via sums.
  # - Prior on alpha: Uniform(0, A). Inside (0, A) its log prior is constant (0),
  #   outside the interval, log prior is -Inf.
  # - Prior on beta: Gamma(gamma0, lambda0) with 'rate' parameterization.
  # - Proposals: independent Gaussian random walks on alpha and beta, consistent
  #   with the reference style. We clip invalid proposals by reverting to current
  #   values, mirroring the reference behavior for sig2.
  #
  calpha <- start[1]; cbeta <- start[2]
  gamma0 <- priorpars[1]; lambda0 <- priorpars[2]; A <- priorpars[3]
  valpha <- jumpvars[1]; vbeta <- jumpvars[2]
  
  alphas <- NULL; betas <- NULL; mus <- NULL
  acceptind <- 0
  cnt <- 0
  
  # Precompute sufficient statistics for the Gamma likelihood
  n <- length(dat)
  sumlogy <- sum(log(dat))
  sumy <- sum(dat)
  
  repeat{
    cnt <- cnt + 1
    alphastar <- proposealpha(calpha, valpha, A)
    betastar  <- proposebeta(cbeta, vbeta)
    
    # log-likelihood (current and proposed)
    #   log f(alpha, beta | y) = n * (alpha * log(beta) - log(Gamma(alpha))) +
    #                            (alpha - 1) * sum(log(y)) - beta * sum(y)
    lfcur  <- n * (calpha * log(cbeta) - lgamma(calpha)) + (calpha - 1) * sumlogy - cbeta * sumy
    lfstar <- n * (alphastar * log(betastar) - lgamma(alphastar)) + (alphastar - 1) * sumlogy - betastar * sumy
    
    # log-prior for alpha: Uniform(0, A)
    #   log pi(alpha) = 0 for alpha in (0, A), -Inf otherwise
    lpi_alpha_cur  <- if(calpha > 0 && calpha < A) 0 else -Inf
    lpi_alpha_star <- if(alphastar > 0 && alphastar < A) 0 else -Inf
    
    # log-prior for beta: Gamma(gamma0, lambda0), rate parameterization
    #   log pi(beta) = gamma0 * log(lambda0) - log(Gamma(gamma0))
    #                  + (gamma0 - 1) * log(beta) - lambda0 * beta
    lpi_beta_cur  <- gamma0 * log(lambda0) - lgamma(gamma0) + (gamma0 - 1) * log(cbeta)  - lambda0 * cbeta
    lpi_beta_star <- gamma0 * log(lambda0) - lgamma(gamma0) + (gamma0 - 1) * log(betastar) - lambda0 * betastar
    
    lpicur  <- lpi_alpha_cur + lpi_beta_cur
    lpistar <- lpi_alpha_star + lpi_beta_star
    
    # Metropolis acceptance (symmetric random-walk proposals)
    astar <- min(exp((lfstar + lpistar) - (lfcur + lpicur)), 1)
    ustar <- runif(1, 0, 1)
    
    newalpha <- calpha; newbeta <- cbeta
    if(ustar <= astar){
      newalpha <- alphastar; newbeta <- betastar
      acceptind <- acceptind + 1
    }
    
    if(cnt > B){
      alphas <- c(alphas, newalpha)
      betas  <- c(betas,  newbeta)
      mus    <- c(mus,    newalpha / newbeta)  # Posterior samples of mu = alpha / beta
    }
    
    calpha <- newalpha; cbeta <- newbeta
    if(cnt == (B + M)) break
  }
  
  cat("acceptprob:", acceptind / M, fill = TRUE)
  res <- data.frame(alpha = alphas, beta = betas, mu = mus)
  attr(res, "acceptprob") <- acceptind / M
  return(res)
}
#----------------------------------------------------------------
proposealpha <- function(calpha, valpha, A){
  # propose jump from random walk for alpha (shape), enforce support (0, A)
  # Reference-style: if invalid, revert to current (like proposesig2 in the ref)
  z <- rnorm(1, 0, sqrt(valpha))
  alphastar <- calpha + z
  if(alphastar <= 0 || alphastar >= A) alphastar <- calpha
  return(alphastar)
}
#----------------------------------------------------------------
proposebeta <- function(cbeta, vbeta){
  # propose jump from random walk for beta (rate), enforce positivity
  # Reference-style: if invalid, revert to current
  z <- rnorm(1, 0, sqrt(vbeta))
  betastar <- cbeta + z
  if(betastar <= 0) betastar <- cbeta
  return(betastar)
}
#----------------------------------------------------------------

### Gibbs 
gibbsforgamma <- function(dat, start, priorpars, B, M, valpha){
  # Gibbs sampler for one-sample Gamma(shape = alpha, rate = beta) model
  # with alpha ~ Uniform(0, A), beta ~ Gamma(gamma0, lambda0)
  #
  # dat        : vector of observed positive data (y_i > 0)
  # start      : c(alpha0, beta0)  -- starting values
  # priorpars  : c(gamma0, lambda0, A)
  # B          : burn-in iterations
  # M          : number of kept Monte Carlo draws
  # valpha     : proposal variance for MH step on alpha (Metropolis-within-Gibbs)
  #
  # Notes on full conditionals and conjugacy:
  # - Conditional for beta | alpha, y is Gamma(gamma0 + n*alpha, lambda0 + sum(y))
  #   (shape/rate parametrization) -- this is conjugate, so we can sample beta directly.
  # - Conditional for alpha | beta, y is NOT standard:
  #       p(alpha | beta, y) proportional to [beta^(n*alpha) / Gamma(alpha)^n] *
  #       (prod(y_i))^(alpha - 1) * I(0 < alpha < A)
  #   We use a random-walk MH step for alpha inside the Gibbs loop
  #   (Metropolis-within-Gibbs), mirroring the reference Gibbs code structure.
  #
  calpha <- start[1]; cbeta <- start[2]
  gamma0 <- priorpars[1]; lambda0 <- priorpars[2]; A <- priorpars[3]
  
  alphas <- NULL; betas <- NULL; mus <- NULL
  cnt <- 0
  accept_alpha <- 0
  
  # Precompute sufficient statistics
  n <- length(dat)
  sumlogy <- sum(log(dat))
  sumy <- sum(dat)
  
  repeat{
    cnt <- cnt + 1
    
    # 1) Sample beta | alpha, y  (conjugate Gamma)
    #    shape = gamma0 + n*alpha ; rate = lambda0 + sum(y)
    newbeta <- rgamma(1, shape = gamma0 + n * calpha, rate = lambda0 + sumy)
    
    # 2) Sample alpha | beta, y  (Metropolis step within Gibbs)
    #    target log-density up to constant:
    #       log p(alpha | beta, y) = n*alpha*log(beta) - n*log(Gamma(alpha))
    #                                + (alpha - 1)*sum(log(y)), for 0<alpha<A
    astep <- sampalpha_mh(calpha, newbeta, valpha, sumlogy, n, A)
    newalpha <- astep$alpha
    accept_alpha <- accept_alpha + astep$acc
    
    if(cnt > B){
      alphas <- c(alphas, newalpha)
      betas  <- c(betas,  newbeta)
      mus    <- c(mus,    newalpha / newbeta)
    }
    calpha <- newalpha; cbeta <- newbeta
    
    if(cnt == (B + M)) break
  }
  
  cat("alpha_acceptprob (within Gibbs):", accept_alpha / M, fill = TRUE)
  res <- data.frame(alpha = alphas, beta = betas, mu = mus)
  return(res)
}
#----------------------------------------------------------------
sampalpha_mh <- function(calpha, beta, valpha, sumlogy, n, A){
  # One-step random-walk MH update for alpha (shape) given beta and y.
  # Returns a list(alpha = ..., acc = 0/1)
  #
  # target log-density (up to constant in alpha):
  #   log f(alpha | beta, y) = n*alpha*log(beta) - n*log(Gamma(alpha))
  #                            + (alpha - 1)*sumlogy
  # with support 0 < alpha < A; outside support, log-density = -Inf.
  #
  z <- rnorm(1, 0, sqrt(valpha))
  alphastar <- calpha + z
  if(alphastar <= 0 || alphastar >= A){
    # As in reference style, invalid proposal -> revert (equivalent to reject)
    return(list(alpha = calpha, acc = 0))
  }
  
  # log target at current and proposed
  lfcur  <- n * calpha   * log(beta) - n * lgamma(calpha)   + (calpha   - 1) * sumlogy
  lfstar <- n * alphastar * log(beta) - n * lgamma(alphastar) + (alphastar - 1) * sumlogy
  
  a <- min(exp(lfstar - lfcur), 1)
  u <- runif(1, 0, 1)
  if(u <= a) return(list(alpha = alphastar, acc = 1))
  return(list(alpha = calpha, acc = 0))
}
#----------------------------------------------------------------