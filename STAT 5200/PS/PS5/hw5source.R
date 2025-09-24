nll_gamma <- function(par, y) {
  a <- exp(par[1]); b <- exp(par[2])      # work on log-scale for positivity
  -sum( a*log(b) - lgamma(a) + (a-1)*log(y) - b*y )
}

fit_by_optim <- function(y) {
  y <- y[y > 0 & is.finite(y)]
  # MoM start on log-scale
  ybar <- mean(y); v <- mean((y-ybar)^2)
  a0 <- if (v > 0) (ybar^2)/v else 1
  b0 <- a0 / ybar
  o <- optim(par = log(c(a0, b0)), fn = nll_gamma, y = y, hessian = TRUE, method = "BFGS")
  ah <- exp(o$par[0+1]); bh <- exp(o$par[1+1])
  J  <- diag(c(ah, bh))                   # Jacobian of exp at the optimum
  cov <- J %*% solve(o$hessian) %*% J     # delta method to natural scale
  se  <- sqrt(diag(cov))
  list(alpha_hat = ah, beta_hat = bh, cov = cov, se_alpha = se[1], se_beta = se[2])
}

## Solve alpha from  log(alpha) - digamma(alpha) = log(mean(y)) - mean(log(y))
solve_alpha_gamma <- function(y, tol = 1e-10, maxit = 100) {
  y <- as.numeric(y)
  y <- y[is.finite(y) & y > 0]
  n <- length(y)
  if (n == 0L) stop("No positive finite observations provided.")
  ybar   <- mean(y)
  clog   <- log(ybar) - mean(log(y))
  
  ## Method-of-moments start: alpha0 = ybar^2 / var
  v <- mean((y - ybar)^2)
  a <- if (v > 0) (ybar^2) / v else 1
  a <- max(a, 1e-6)
  
  for (iter in 1:maxit) {
    f <- log(a) - digamma(a) - clog
    g <- 1/a - trigamma(a)      # derivative of f
    step <- f / g
    a_new <- a - step
    if (!is.finite(a_new) || a_new <= 0) a_new <- a / 2
    if (abs(a_new - a) < tol * (1 + abs(a))) return(a_new)
    a <- a_new
  }
  warning("Alpha solver reached maxit without tight convergence; returning last iterate.")
  a
}

## MLE + Wald SEs/CIs for Gamma(shape \alpha, rate \beta)
gamma_mle_wald <- function(y, conf = 0.95) {
  y <- as.numeric(y)
  y <- y[is.finite(y) & y > 0]
  n <- length(y)
  if (n == 0L) stop("No positive finite observations provided.")
  
  ybar     <- mean(y)
  sumy     <- sum(y)
  sumlogy  <- sum(log(y))
  
  ## MLEs
  alpha_hat <- solve_alpha_gamma(y)
  beta_hat  <- alpha_hat / ybar
  
  ## Observed information (shape–rate):
  ## I11 = n * trigamma(alpha), I22 = n * alpha / beta^2, I12 = - n / beta
  I11 <- n * trigamma(alpha_hat)
  I22 <- n * alpha_hat / (beta_hat^2)
  I12 <- - n / beta_hat
  I   <- rbind(c(I11, I12), c(I12, I22))
  
  cov2 <- solve(I)                     # inverse information
  se_a <- sqrt(cov2[1,1])
  se_b <- sqrt(cov2[2,2])
  
  ## Wald CIs on natural scale
  z <- qnorm(0.5 + conf/2)
  ci_a <- c(alpha_hat - z*se_a, alpha_hat + z*se_a)
  ci_b <- c(beta_hat  - z*se_b, beta_hat  + z*se_b)
  
  ## Log-scale Wald CIs (delta method) for positivity
  se_log_a <- se_a / alpha_hat
  se_log_b <- se_b / beta_hat
  ci_log_a <- exp(log(alpha_hat) + c(-1,1)*z*se_log_a)
  ci_log_b <- exp(log(beta_hat)  + c(-1,1)*z*se_log_b)
  
  list(
    n = n,
    ybar = ybar,
    alpha_hat = alpha_hat,
    beta_hat  = beta_hat,
    se_alpha = se_a,
    se_beta  = se_b,
    Wald95_alpha = ci_a,
    Wald95_beta  = ci_b,
    logWald95_alpha = ci_log_a,
    logWald95_beta  = ci_log_b,
    cov = cov2,
    suff_stats = list(sum_y = sumy, sum_log_y = sumlogy)
  )
}

## Tidy summary data frame for your report/table (no “screen dumps” needed)
to_row <- function(name, r) {
  data.frame(
    group = name,
    n = r$n,
    ybar = r$ybar,
    alpha_hat = r$alpha_hat,
    se_alpha  = r$se_alpha,
    Wald95_alpha_L = r$Wald95_alpha[1],
    Wald95_alpha_U = r$Wald95_alpha[2],
    beta_hat  = r$beta_hat,
    se_beta   = r$se_beta,
    Wald95_beta_L  = r$Wald95_beta[1],
    Wald95_beta_U  = r$Wald95_beta[2],
    logWald95_alpha_L = r$logWald95_alpha[1],
    logWald95_alpha_U = r$logWald95_alpha[2],
    logWald95_beta_L  = r$logWald95_beta[1],
    logWald95_beta_U  = r$logWald95_beta[2]
  )
}

# Helper to solve for alpha (shape) in Gamma(shape \alpha, rate \beta)
solve_alpha_gamma <- function(y, tol = 1e-10, maxit = 100) {
  y <- as.numeric(y); y <- y[is.finite(y) & y > 0]
  ybar <- mean(y); clog <- log(ybar) - mean(log(y))
  v <- mean((y - ybar)^2)
  a <- if (v > 0) (ybar^2)/v else 1
  a <- max(a, 1e-6)
  for (it in 1:maxit) {
    f <- log(a) - digamma(a) - clog
    g <- 1/a - trigamma(a)
    a_new <- a - f/g
    if (!is.finite(a_new) || a_new <= 0) a_new <- a/2
    if (abs(a_new - a) < tol * (1 + abs(a))) return(a_new)
    a <- a_new
  }
  warning("alpha solver hit maxit"); a
}

mle_gamma_sr <- function(y) {
  y <- y[y > 0 & is.finite(y)]
  a <- solve_alpha_gamma(y)
  b <- a / mean(y)
  list(alpha = a, beta = b)
}

loglik_gamma_sr <- function(y, a, b) {
  y <- y[y > 0 & is.finite(y)]
  n <- length(y)
  n * (a * log(b) - lgamma(a)) + (a - 1) * sum(log(y)) - b * sum(y)
}

fit_gamma_sr <- function(y) {
  y <- as.numeric(y); y <- y[is.finite(y) & y > 0]
  if (!length(y)) stop("No positive finite observations.")
  n <- length(y)
  ybar <- mean(y)
  
  a_hat <- solve_alpha_gamma(y)
  b_hat <- a_hat / ybar
  
  ## Observed information at MLE:
  ## I11 = n * trigamma(alpha),  I22 = n * alpha / beta^2,  I12 = - n / beta
  I11 <- n * trigamma(a_hat)
  I22 <- n * a_hat / (b_hat^2)
  I12 <- - n / b_hat
  I   <- rbind(c(I11, I12), c(I12, I22))
  
  cov_ab <- solve(I)  # Var[alpha, beta]
  list(alpha = a_hat, beta = b_hat, cov = cov_ab, n = n, ybar = ybar)
}

## Delta-method mean = alpha/beta, gradient = (1/beta, -alpha/beta^2)
mean_from_fit <- function(fit, conf = 0.95) {
  a <- fit$alpha; b <- fit$beta; S <- fit$cov
  mu_hat <- a / b
  g <- c(1/b, -a/(b^2))
  var_mu <- as.numeric(t(g) %*% S %*% g)
  se_mu  <- sqrt(var_mu)
  z <- qnorm(0.5 + conf/2)
  ci_mu <- c(mu_hat - z*se_mu, mu_hat + z*se_mu)
  list(mu = mu_hat, se = se_mu, ci = ci_mu)
}

# Helper to summarize a t.test object + group descriptives
summarize_ttest <- function(y1, y2, label, var_equal = FALSE) {
  y1 <- y1[is.finite(y1)]
  y2 <- y2[is.finite(y2)]
  tt <- t.test(y1, y2, var.equal = var_equal)  # Welch by default unless var.equal=TRUE
  out <- data.frame(
    analysis     = label,
    equal_var    = var_equal,
    n1           = length(y1),
    n2           = length(y2),
    mean1        = mean(y1),
    mean2        = mean(y2),
    sd1          = sd(y1),
    sd2          = sd(y2),
    t_stat       = unname(tt$statistic),
    df           = unname(tt$parameter),
    p_value      = tt$p.value,
    ci_lower     = tt$conf.int[1],
    ci_upper     = tt$conf.int[2],
    method       = tt$method
  )
  out
}

# Delta-method for the mode m = (alpha - 1)/beta (valid when alpha > 1)
mode_from_fit <- function(fit, conf = 0.95) {
  a <- fit$alpha; b <- fit$beta; S <- fit$cov
  if (a <= 1) {
    # Boundary case: mode at 0. Delta-method on boundary can misbehave.
    # We report mode=0 and omit a symmetric Wald CI (or could give a one-sided CI).
    return(list(mode = 0, se = NA_real_, ci = c(NA_real_, NA_real_), boundary = TRUE))
  }
  m_hat <- (a - 1) / b
  g <- c(1/b, -(a - 1)/(b^2))
  var_m <- as.numeric(t(g) %*% S %*% g)
  se_m  <- sqrt(var_m)
  z <- qnorm(0.5 + conf/2)
  ci_m <- c(m_hat - z*se_m, m_hat + z*se_m)
  list(mode = m_hat, se = se_m, ci = ci_m, boundary = FALSE)
}