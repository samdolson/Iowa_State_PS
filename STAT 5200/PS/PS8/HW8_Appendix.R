# Note:
# The following functions were written to perform the primary computations used in this assignment
# Posterior predictive p-values and credible intervals are computed within the body of the analysis 
# (I didn't take a purely functional approach) 
# 
# To (attempt to) ease your grading, you can keyterm search the following for those calculations/methods 
# 5-number Summaries, Credible Interval Calculations: "fivenum", "ci_alpha", "ci_beta", etc. 
# Predictive p-values: "ppp", "yrep_stats" 

library(knitr)
library(kableExtra)

metropforgamma <- function(dat, start, priorpars, jumpvars, B, M){
  # Metropolis for a one-sample Gamma(shape = alpha, rate = beta) model
  #   with joint prior = product of marginals 
  # alpha ~ Uniform(0, A)
  # beta ~ Gamma(gamma0, lambda0)
  #
  # Inputs: 
  # dat : vector of observed positive data (y_i > 0)
  # start : starting values for (alpha, beta) c(alpha0, beta0)
  # priorpars: prior parameters c(gamma0, lambda0, A)
  #   gamma0, lambda0 are shape/rate of prior on beta
  #   A is the upper bound for alpha's Uniform(0, A) prior
  # jumpvars : proposal variances for random-walk jumps c(valpha, vbeta)
  # B : burn-in iterations
  # M : number of kept draws
  #
  # Additional Notes:
  # Likelihood: Y_i ~ Gamma(alpha, beta)
  #   f(y | alpha, beta) = beta^alpha / Gamma(alpha) * y^(alpha-1) * exp(-beta*y)
  # Prior on alpha: Uniform(0, A) 
  # Prior on beta: Gamma(gamma0, lambda0) also using 'rate' parameterization
  # Proposals: independent Gaussian random walks on alpha and beta, consistent
  #   Adjust invalid proposals (negative parameter values) by reverting to current values
  #   I.e., use same method as the Normal example from Kaiser for sig2
  
  # Starting inputs 
  calpha <- start[1]
  cbeta <- start[2]
  gamma0 <- priorpars[1]
  lambda0 <- priorpars[2]
  A <- priorpars[3]
  valpha <- jumpvars[1]
  vbeta <- jumpvars[2]
  
  # Initialize 
  alphas <- NULL
  betas <- NULL
  mus <- NULL
  
  acceptind <- 0
  cnt <- 0
  
  # Calculate sufficient statistics for the Gamma likelihood
  # Based on initial data input 
  n <- length(dat)
  sumlogy <- sum(log(dat))
  sumy <- sum(dat)
  
  repeat{
    cnt <- cnt + 1
    alphastar <- proposealpha(calpha, valpha, A)
    betastar <- proposebeta(cbeta, vbeta)
    
    # log-likelihood 
    # log f(alpha, beta | y) = n * (alpha * log(beta) - log(Gamma(alpha))) + (alpha - 1) * sum(log(y)) - beta * sum(y)
    lfcur <- n * (calpha * log(cbeta) - lgamma(calpha)) + 
      (calpha - 1) * sumlogy - cbeta * sumy
    lfstar <- n * (alphastar * log(betastar) - lgamma(alphastar)) + 
      (alphastar - 1) * sumlogy - betastar * sumy
    
    # log-prior for alpha: Uniform(0, A)
    #   log pi(alpha) = 0 for alpha in (0, A), -Inf otherwise
    lpi_alpha_cur <- if(calpha > 0 && calpha < A) 0 else -Inf
    lpi_alpha_star <- if(alphastar > 0 && alphastar < A) 0 else -Inf
    
    # log-prior for beta: Gamma(gamma0, lambda0), rate parameterization
    # log pi(beta) = gamma0 * log(lambda0) - log(Gamma(gamma0)) + (gamma0 - 1) * log(beta) - lambda0 * beta
    lpi_beta_cur <- gamma0 * log(lambda0) - lgamma(gamma0) + 
      (gamma0 - 1) * log(cbeta)  - lambda0 * cbeta
    lpi_beta_star <- gamma0 * log(lambda0) - lgamma(gamma0) + 
      (gamma0 - 1) * log(betastar) - lambda0 * betastar
    
    lpicur <- lpi_alpha_cur + lpi_beta_cur
    lpistar <- lpi_alpha_star + lpi_beta_star
    
    # Metropolis acceptance calculations
    astar <- min(exp((lfstar + lpistar) - (lfcur + lpicur)), 1)
    ustar <- runif(1, 0, 1)
    
    newalpha <- calpha
    newbeta <- cbeta
    if(ustar <= astar){
      # update
      newalpha <- alphastar
      newbeta <- betastar
      acceptind <- acceptind + 1
    }
    
    if(cnt > B){
      alphas <- c(alphas, newalpha)
      betas <- c(betas, newbeta)
      # Posterior samples of mu 
      # Generate mu based on the new alpha AND beta
      mus <- c(mus, newalpha / newbeta)  
    }
    
    calpha <- newalpha
    cbeta <- newbeta
    if(cnt == (B + M)) break
  }
  
  cat("acceptprob:", acceptind / M, fill = TRUE)
  res <- data.frame(alpha = alphas, beta = betas, mu = mus)
  attr(res, "acceptprob") <- acceptind / M
  return(res)
}

proposealpha <- function(calpha, valpha, A){
  # propose jump from random walk for alpha (shape)
  # check support in (0, A)
  # if invalid, revert to current
  z <- rnorm(1, 0, sqrt(valpha))
  alphastar <- calpha + z
  if(alphastar <= 0 || alphastar >= A) alphastar <- calpha
  return(alphastar)
}

proposebeta <- function(cbeta, vbeta){
  # propose jump from random walk for beta (rate)
  # check positivity
  # if invalid, revert to current
  z <- rnorm(1, 0, sqrt(vbeta))
  betastar <- cbeta + z
  if(betastar <= 0) betastar <- cbeta
  return(betastar)
}

gibbsforgamma <- function(dat, start, priorpars, B, M, valpha){
  # Gibbs sampler for one-sample Gamma(shape = alpha, rate = beta) model
  #   with alpha ~ Uniform(0, A)
  #   beta ~ Gamma(gamma0, lambda0)
  # 
  # Inputs: 
  #   (similar to MH method)
  # dat : vector of observed positive data (y_i > 0)
  # start : c(alpha0, beta0)
  # priorpars : c(gamma0, lambda0, A)
  # B : burn-in 
  # M : number of kept draws
  # valpha : proposal variance for MH step on alpha (bc Metropolis-within-Gibbs)
  #
  # Notes on full conditionals and conjugacy:
  # Conditional for beta | alpha, y is Gamma(gamma0 + n*alpha, lambda0 + sum(y))
  #   (shape/rate parametrization) this is conjugate, so we can sample beta directly
  # Conditional for alpha | beta, y (the tricky part):
  #       p(alpha | beta, y) proportional to [beta^(n*alpha) / Gamma(alpha)^n] *
  #       (prod(y_i))^(alpha - 1) * I(0 < alpha < A)
  # Use a random-walk MH for alpha inside the Gibbs
  
  # Starting inputs 
  calpha <- start[1]
  cbeta <- start[2]
  gamma0 <- priorpars[1]
  lambda0 <- priorpars[2]
  A <- priorpars[3]
  
  # Initialize 
  alphas <- NULL
  betas <- NULL
  mus <- NULL
  
  cnt <- 0
  accept_alpha <- 0
  
  # Same sufficient statistics as before
  n <- length(dat)
  sumlogy <- sum(log(dat))
  sumy <- sum(dat)
  
  repeat{
    cnt <- cnt + 1
    
    # 1. Sample beta | alpha, y  (conjugate Gamma)
    #   where 
    #   shape = gamma0 + n*alpha
    #   rate = lambda0 + sum(y)
    newbeta <- rgamma(1, shape = gamma0 + n * calpha, rate = lambda0 + sumy)
    
    # 2. Sample alpha | beta, y  (Metropolis within Gibbs)
    #   target log-density up to constant:
    #   log p(alpha | beta, y) = n*alpha*log(beta) - n*log(Gamma(alpha)) + (alpha - 1)*sum(log(y))
    #   for 0<alpha<A
    astep <- sampalpha_mh(calpha, newbeta, valpha, sumlogy, n, A)
    newalpha <- astep$alpha
    accept_alpha <- accept_alpha + astep$acc
    
    if(cnt > B){
      alphas <- c(alphas, newalpha)
      betas <- c(betas,  newbeta)
      mus <- c(mus,    newalpha / newbeta)
    }
    calpha <- newalpha
    cbeta <- newbeta
    
    if(cnt == (B + M)) break
  }
  
  cat("alpha_acceptprob (within Gibbs):", accept_alpha / M, fill = TRUE)
  res <- data.frame(alpha = alphas, beta = betas, mu = mus)
  return(res)
}

sampalpha_mh <- function(calpha, beta, valpha, sumlogy, n, A){
  # MH for alpha given beta and y.
  # 
  # target log-density (up to constant in alpha):
  #   log f(alpha | beta, y) = n*alpha*log(beta) - n*log(Gamma(alpha)) + (alpha - 1)*sumlogy
  #   with support 0 < alpha < A
  z <- rnorm(1, 0, sqrt(valpha))
  alphastar <- calpha + z
  if(alphastar <= 0 || alphastar >= A){
    # If proposal is invalid
    # revert to prior value
    return(list(alpha = calpha, acc = 0))
  }
  
  # log target at current and proposed
  lfcur <- n * calpha   * log(beta) - n * lgamma(calpha)   + 
    (calpha   - 1) * sumlogy
  lfstar <- n * alphastar * log(beta) - n * lgamma(alphastar) + 
    (alphastar - 1) * sumlogy
  
  a <- min(exp(lfstar - lfcur), 1)
  u <- runif(1, 0, 1)
  if(u <= a) return(list(alpha = alphastar, acc = 1))
  return(list(alpha = calpha, acc = 0))
}

# Done in conjunction w/Burn-In evaluation 
runmean <- function(x) cumsum(x) / seq_along(x)

# Done to adjust table formatting 
esc <- function(x) gsub("%", "\\%", x, fixed = TRUE)

one_method_tbl <- function(vec_list, method = NULL) {
  out <- do.call(
    rbind,
    lapply(vec_list, function(x) {
      fn <- fivenum(x)
      ci <- quantile(x, c(0.025, 0.975), names = FALSE)
      c(Min = fn[1], Q1 = fn[2], Median = fn[3], Q3 = fn[4], Max = fn[5],
        "CI 2.5" = ci[1], "CI 97.5" = ci[2])
    })
  )
  df <- data.frame(Parameter = names(vec_list), round(out, 4),
                   row.names = NULL)
  if (!is.null(method)) df <- cbind(Method = method, df)
  df
}

# Posterior Predictive p-value statistics 
#   these can include NA handling, but non-issue for this problem
Q1 <- function(z) as.numeric(quantile(z, 0.75))
Q2 <- function(z) diff(range(z))

# For Gibbs Tuning 
alpha_accept_from_chain <- function(alpha_chain) {
  if (length(alpha_chain) < 2) return(NA)
  mean(alpha_chain[-1] != alpha_chain[-length(alpha_chain)])
}

# For Question 6 plot comparisons
make_breaks <- function(x1, x2, width = 0.1) seq(from = floor(x1), to = ceiling(x2), by = width)

# Actual Run of Code 
# If you're just interested in functions, you can stop here, well, kinda 
# gammaDat <- read.table("C:/Users/samue/OneDrive/Desktop/Iowa_State_PS/STAT 5200/PS/PS8/gammadat_bayes.txt", header = T)
gammaDat <- read.table(".../gammadat_bayes.txt", header = T)
set.seed(43)
y <- as.numeric(gammaDat$y)
m <- mean(y)
v <- var(y)

# MH 
alpha0 <- if (v > 0) (m^2 / v) else 1.0
alpha0 <- max(0.10, min(alpha0, 19.9))
beta0 <- max(1e-6, alpha0 / max(m, 1e-6))
start1 <- c(alpha0, beta0)

valpha <- (0.20 * alpha0)^2
vbeta <- (0.20 * beta0)^2
jumpvars1 <- c(valpha, vbeta)

set.seed(43)
j_diag <- metropforgamma(
  dat = y,
  start = start1,
  priorpars = c(gamma0 = 0.5, lambda0 = 0.1, A = 20),
  jumpvars = jumpvars1,
  B = 0,
  M = 50000
)

B_line <- 5000L   

par(mfrow = c(1,2))
plot(j_diag$alpha, type="l", main="alpha trace (B=0)", xlab="iter")
abline(v = B_line, col=2, lty=2, lwd=2)
abline(v = B_line*2, col=2, lty=2, lwd=2)
lines(runmean(j_diag$alpha), col = 4)

plot(j_diag$beta,  type="l", main="beta trace (B=0)",  xlab="iter")
abline(v = B_line, col=2, lty=2, lwd=2)
abline(v = B_line*2, col=2, lty=2, lwd=2)
lines(runmean(j_diag$beta),  col = 4)

par(mfrow = c(1,2))
acf(j_diag$alpha, lag.max = 5000, main="ACF alpha (full)")
acf(j_diag$beta,  lag.max = 5000, main="ACF beta (full)")

scales <- c(0.5, 0.75, 1.0, 1.25, 1.5)
set.seed(43)
tune_tab <- do.call(rbind, lapply(scales, function(s){
  j <- metropforgamma(
    dat = y, start = start1,
    priorpars = c(0.5, 0.1, 20),
    jumpvars = jumpvars1 * s^2,
    B = 10000, M = 50000
  )
  data.frame(scale = s,
             valpha = (jumpvars1[1] * s^2),
             vbeta = (jumpvars1[2] * s^2),
             accept = attr(j, "acceptprob"))
}))

# Final MH 
set.seed(43)
final <- metropforgamma(
  dat = y,
  start = start1,
  priorpars = c(0.5, 0.1, 20),
  # from tuning
  jumpvars = jumpvars1 * 1.25^2,   
  B = 10000,                       
  M = 50000
)

# Gather relevant statistics from final
acc_final <- attr(final, "acceptprob")
alpha <- final$alpha
beta <- final$beta
mu <- final$mu
fivenum_alpha <- fivenum(alpha)
ci_alpha <- quantile(alpha, c(0.025, 0.975))
fivenum_beta <- fivenum(beta)
ci_beta <- quantile(beta,  c(0.025, 0.975))
fivenum_mu <- fivenum(mu)
ci_mu <- quantile(mu,    c(0.025, 0.975))
corr_ab <- cor(alpha, beta, use = "complete.obs")

alphaMH <- alpha
betaMH <- beta
muMH <- mu
ci_alphaMH <- ci_alpha
ci_betaMH <- ci_beta
ci_muMH <- ci_mu

mh_tbl <- one_method_tbl(list(alpha = alphaMH, beta = betaMH, mu = muMH), method = "MH")

if (exists("alphaG") && exists("betaG") && exists("muG")) {
  gibbs_tbl <- one_method_tbl(list(alpha = alphaG, beta = betaG, mu = muG), method = "Gibbs")
  param_tbl <- rbind(mh_tbl, gibbs_tbl)
} else {
  param_tbl <- mh_tbl
}

kable(
  param_tbl, booktabs = TRUE, escape = TRUE,
  caption = esc("Posterior summaries: five-number statistics and 95 central credible intervals")
) |>
  kable_styling(full_width = FALSE, position = "center",
                latex_options = c("hold_position")) |>
  column_spec(1, bold = TRUE) |>
  add_header_above(c(" " = if ("Method" %in% names(param_tbl)) 2 else 1,
                     "Posterior Summary" = ncol(param_tbl) -
                       if ("Method" %in% names(param_tbl)) 2 else 1))

nm <- names(param_tbl)

ci_low <- nm[grepl("^CI(\\.| )2\\.5",  nm)]
ci_high <- nm[grepl("^CI(\\.| )97\\.5", nm)]

cols_ci <- c(if ("Method" %in% nm) "Method", "Parameter", ci_low, ci_high)
ci_only <- param_tbl[, cols_ci, drop = FALSE]

names(ci_only)[names(ci_only) == ci_low] <- "Lower 2.5"
names(ci_only)[names(ci_only) == ci_high] <- "Upper 97.5"

kable(
  ci_only, booktabs = TRUE, escape = TRUE,
  caption = esc("95 central credible intervals")
) |>
  kable_styling(full_width = FALSE, position = "center",
                latex_options = c("hold_position")) |>
  column_spec(1, bold = TRUE)

acc_vec <- as.numeric(acc_final)
acc_names <- names(acc_final)
if (is.null(acc_names)) {
  acc_names <- if (length(acc_vec) == 1) "Overall" else paste0("Component_", seq_along(acc_vec))
}
acc_tbl <- data.frame(Component = acc_names,
                      `Acceptance Rate` = paste0(round(100 * acc_vec, 1), "%"),
                      row.names = NULL)

kable(
  acc_tbl, booktabs = TRUE, escape = TRUE,
  caption = "Acceptance rates"
) |>
  kable_styling(full_width = FALSE, position = "center",
                latex_options = c("hold_position")) |>
  column_spec(1, bold = TRUE)

corr_mat <- matrix(c(1, corr_ab, corr_ab, 1), nrow = 2,
                   dimnames = list(c("alpha","beta"), c("alpha","beta")))
corr_df <- data.frame(Parameter = rownames(corr_mat), round(corr_mat, 4), row.names = NULL)

kable(
  corr_df, booktabs = TRUE, escape = TRUE,
  caption = "Correlation matrix for alpha and beta"
) |>
  kable_styling(full_width = FALSE, position = "center",
                latex_options = c("hold_position")) |>
  column_spec(1, bold = TRUE)

par(mfrow = c(1,3))
hist(alphaMH, breaks="FD", main=expression(paste("Posterior of ", alpha)), xlab=expression(alpha))
abline(v=ci_alpha, lty=2)
hist(betaMH, breaks="FD", main=expression(paste("Posterior of ", beta)),  xlab=expression(beta))
abline(v=ci_beta,  lty=2)
hist(muMH, breaks="FD", main=expression(paste("Posterior of ", mu==alpha/beta)), xlab=expression(mu))
abline(v=ci_mu, lty=2)

# Quartile obs
Q_obs1 <- Q1(y)
# Range obs
Q_obs2 <- Q2(y)

# 10,000 draws 
ndraws <- 10000L
n <- length(y)

set.seed(43)
# randomly take the 10,000
idx <- sample(seq_along(alpha), ndraws)

yrep_stats <- matrix(NA, nrow = 2, ncol = ndraws)

# for each of the 10,000, create datasets of original size of obs (50)
for (s in seq_len(ndraws)) {
  yrep <- rgamma(n, shape = alpha[idx[s]], rate = beta[idx[s]])
  yrep_stats[1, s] <- Q1(yrep)  
  yrep_stats[2, s] <- Q2(yrep)  
}

# upper-tail for 75th percentile
ppp_Q75 <- mean(yrep_stats[1, ] >= Q_obs1)  
# upper-tail for range
ppp_range <- mean(yrep_stats[2, ] >= Q_obs2)  

ppp_tbl <- data.frame(
  Statistic = c("75th percentile", "Range"),
  Observed = round(c(Q_obs1, Q_obs2), 4),
  PPP = round(c(ppp_Q75, ppp_range), 4),
  check.names = FALSE
)

kable(
  ppp_tbl, booktabs = TRUE,
  caption = "Posterior predictive p-values (upper) 10,000 datasets"
) |>
  kable_styling(full_width = FALSE, position = "center",
                latex_options = c("hold_position")) |>
  column_spec(1, bold = TRUE)

# Gibbs
A <- 20 
alpha0 <- max(0.10, min(ifelse(v > 0, m^2 / v, 1.0), A - 0.1))
beta0 <- max(1e-6, alpha0 / max(m, 1e-6))
start <- c(alpha0, beta0)
priorpars <- c(gamma0 = 0.5, lambda0 = 0.1, A = A)

# Metropolis-within-Gibbs
# only alpha uses MH
valpha <- (0.15 * alpha0)^2

set.seed(43)
diag_all <- gibbsforgamma(
  dat = y, start = start, priorpars = priorpars,
  B = 0, M = 50000, valpha = valpha
)

alpha_all <- diag_all$alpha
beta_all <- diag_all$beta
mu_all <- diag_all$mu

B_line <- 10000L

par(mfrow = c(1,2))
plot(alpha_all, type="l", main="alpha trace (B=0)", xlab="iter", ylab=expression(alpha))
abline(v=B_line, col=2, lty=2, lwd=2)
abline(v=B_line*2, col=2, lty=2, lwd=2)
lines(runmean(alpha_all), col=4)

plot(beta_all,  type="l", main="beta trace (B=0)",  xlab="iter", ylab=expression(beta))
abline(v=B_line, col=2, lty=2, lwd=2)
abline(v=B_line*2, col=2, lty=2, lwd=2)
lines(runmean(beta_all),  col=4)

par(mfrow = c(1,2))
acf(alpha_all, lag.max=10000, main="ACF alpha (full)")
acf(beta_all,  lag.max=10000, main="ACF beta (full)")

scales <- c(0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
set.seed(43)
tune_tab <- do.call(rbind, lapply(scales, function(s){
  j <- gibbsforgamma(
    dat = y, start = start, priorpars = priorpars,
    B = 20000, M = 50000,
    valpha = valpha * s^2
  )
  data.frame(
    scale = s,
    valpha = valpha * s^2,
    acc_alpha = if (!is.null(j$alpha_accept_overall))
      j$alpha_accept_overall
    else
      alpha_accept_from_chain(j$alpha)
  )
}))

# Final Gibbs 
B_line <- 20000L 

set.seed(43)
final <- gibbsforgamma(
  dat = y, start = start, priorpars = priorpars,
  B = B_line, M = 50000,
  # tuning parameter
  valpha = valpha * 2^2   
)

# Visuals and analysis for Gibbs 
acc_alpha_final <- if (!is.null(final$alpha_accept_overall)) final$alpha_accept_overall else NA

alpha <- final$alpha
beta <- final$beta
mu <- final$mu 

fivenum_alpha <- fivenum(alpha)
fivenum_beta <- fivenum(beta)
ci_alpha <- quantile(alpha, c(0.025, 0.975), names = FALSE)
ci_beta <- quantile(beta,  c(0.025, 0.975), names = FALSE)
ci_mu <- quantile(mu, c(0.025, 0.975), names = FALSE) 

corr_ab <- cor(alpha, beta, use = "complete.obs")

# "Credible Interval" 
#   Is based on empirical distribution of the parameter(s)
summ_tbl <- data.frame(
  Parameter = c("alpha","beta"),
  Min = c(fivenum_alpha[1], fivenum_beta[1]),
  Q1 = c(fivenum_alpha[2], fivenum_beta[2]),
  Median = c(fivenum_alpha[3], fivenum_beta[3]),
  Q3 = c(fivenum_alpha[4], fivenum_beta[4]),
  Max = c(fivenum_alpha[5], fivenum_beta[5]),
  `CI 2.5%` = c(ci_alpha[1], ci_beta[1]),
  `CI 97.5%` = c(ci_alpha[2], ci_beta[2]),
  check.names = FALSE
)
summ_tbl[,-1] <- lapply(summ_tbl[,-1], function(z) round(as.numeric(z), 4))

kable(
  summ_tbl, booktabs = TRUE, escape = TRUE,
  caption = esc("Posterior summaries for alpha and beta: five-number statistics and 95 central credible intervals")
) |>
  kable_styling(full_width = FALSE, position = "center",
                latex_options = c("hold_position")) |>
  column_spec(1, bold = TRUE) |>
  add_header_above(c(" " = 1, "Posterior Summary" = ncol(summ_tbl) - 1))

corr_mat <- matrix(c(1, corr_ab, corr_ab, 1), nrow = 2,
                   dimnames = list(c("alpha","beta"), c("alpha","beta")))
corr_df <- data.frame(Parameter = rownames(corr_mat), round(corr_mat, 4), row.names = NULL)

kable(
  corr_df, booktabs = TRUE, escape = TRUE,
  caption = "Correlation matrix for alpha and beta (from Markov chain)"
) |>
  kable_styling(full_width = FALSE, position = "center",
                latex_options = c("hold_position")) |>
  column_spec(1, bold = TRUE)

alphaG <- alpha
ci_alphaG <- ci_alpha
betaG <- beta 
ci_betaG <- ci_beta
muG <- mu 
ci_muG <- ci_mu

par(mfrow = c(1,3))
hist(alphaG, breaks="FD", main=expression(paste("Posterior of ", alpha)), xlab=expression(alpha))
abline(v=ci_alpha, lty=2)
hist(betaG,  breaks="FD", main=expression(paste("Posterior of ", beta)),  xlab=expression(beta))
abline(v=ci_beta,  lty=2)
hist(muG,    breaks="FD", main=expression(paste("Posterior of ", mu==alpha/beta)), xlab=expression(mu))
abline(v=ci_mu, lty=2)

gout <- final
n <- length(y)

# Posterior Predictive Stuff 
#   Should be closely mirrored when done for Gibbs (though, maybe slightly different)
# observed
Tobs_q75 <- Q1(y)
Tobs_range <- Q2(y)

# posterior draws
ndraws <- 10000L
K <- nrow(gout)
idx <- if (K >= ndraws) sample.int(K, ndraws) else sample.int(K, ndraws, replace = TRUE)
a_draws <- gout$alpha[idx]
b_draws <- gout$beta[idx]

# posterior predictive datasets
Trep_q75 <- numeric(ndraws)
Trep_range <- numeric(ndraws)
for (k in seq_len(ndraws)) {
  yrep <- rgamma(n, shape = a_draws[k], rate = b_draws[k])
  Trep_q75[k] <- Q1(yrep)
  Trep_range[k] <- Q2(yrep)
}

ppp_q75 <- mean(Trep_q75   >= Tobs_q75)
ppp_range <- mean(Trep_range >= Tobs_range)

ppp_tbl <- data.frame(
  Statistic = c("75th percentile", "Range"),
  Observed = round(c(Tobs_q75, Tobs_range), 4),
  PPP = round(c(ppp_q75, ppp_range), 4),
  check.names = FALSE
)

kable(
  ppp_tbl, booktabs = TRUE,
  caption = "Posterior predictive p-values (upper) 10,000 datasets"
) |>
  kable_styling(full_width = FALSE, position = "center",
                latex_options = c("hold_position")) |>
  column_spec(1, bold = TRUE) |>
  add_header_above(c(" " = 1, "Posterior Predictive Check" = 2))

# Question 6 then was just using some combination of previously created objects and formatting them
alpha_range <- range(c(alphaMH, alphaG), na.rm = TRUE)
beta_range <- range(c(betaMH,  betaG),  na.rm = TRUE)
mu_range <- range(c(muMH,    muG),    na.rm = TRUE)

alpha_breaks <- make_breaks(min(alpha_range), max(alpha_range), width = 0.2)
beta_breaks <- make_breaks(min(beta_range),  max(beta_range),  width = 0.05)
mu_breaks <- make_breaks(min(mu_range),    max(mu_range),    width = 0.1)

alpha_ylim <- c(0, max(
  hist(alphaMH, breaks = alpha_breaks, plot = FALSE)$counts,
  hist(alphaG,  breaks = alpha_breaks, plot = FALSE)$counts
))
beta_ylim <- c(0, max(
  hist(betaMH, breaks = beta_breaks, plot = FALSE)$counts,
  hist(betaG,  breaks = beta_breaks, plot = FALSE)$counts
))
mu_ylim <- c(0, max(
  hist(muMH, breaks = mu_breaks, plot = FALSE)$counts,
  hist(muG,  breaks = mu_breaks, plot = FALSE)$counts
))

par(mfrow = c(1,2))

hist(alphaMH, breaks = alpha_breaks, xlim = alpha_range, ylim = alpha_ylim,
     main = expression(paste("MH - Posterior of ", alpha)),
     xlab = expression(alpha))
abline(v = ci_alphaMH, lty = 2)

hist(alphaG, breaks = alpha_breaks, xlim = alpha_range, ylim = alpha_ylim,
     main = expression(paste("Gibbs - Posterior of ", alpha)),
     xlab = expression(alpha))
abline(v = ci_alphaG, lty = 2)

par(mfrow = c(1,2))

hist(betaMH, breaks = beta_breaks, xlim = beta_range, ylim = beta_ylim,
     main = expression(paste("MH - Posterior of ", beta)),
     xlab = expression(beta))
abline(v = ci_betaMH, lty = 2)

hist(betaG, breaks = beta_breaks, xlim = beta_range, ylim = beta_ylim,
     main = expression(paste("Gibbs - Posterior of ", beta)),
     xlab = expression(beta))
abline(v = ci_betaG, lty = 2)

par(mfrow = c(1,2))

hist(muMH, breaks = mu_breaks, xlim = mu_range, ylim = mu_ylim,
     main = expression(paste("MH - Posterior of ", mu == alpha/beta)),
     xlab = expression(mu))
abline(v = ci_muMH, lty = 2)

hist(muG, breaks = mu_breaks, xlim = mu_range, ylim = mu_ylim,
     main = expression(paste("Gibbs - Posterior of ", mu == alpha/beta)),
     xlab = expression(mu))
abline(v = ci_muG, lty = 2)