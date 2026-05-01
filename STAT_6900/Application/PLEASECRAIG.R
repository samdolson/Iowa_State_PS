# Setup 
library(xts)
library(dplyr)
library(zoo)

set.seed(42)
FREQ  <- 52L;  T_HOR <- 10.0;  DT <- 1/FREQ
K     <- round(T_HOR / DT)
X0    <- 1.0;  Z <- 2.0;  LAM <- 1.0;  M <- 2L

dax            <- readRDS("dax.rds")
dax_clean      <- readRDS("dax_clean.rds")
dax_ret        <- readRDS("dax_ret.rds")
GDAXI          <- readRDS("GDAXI.rds")
r_series       <- readRDS("r_series.rds")
train_ret      <- readRDS("train_ret.rds")
test_ret       <- readRDS("test_ret.rds")
price_list     <- readRDS("price_list.rds")
train_paths_wk <- readRDS("train_paths_wk.rds")
test_paths_wk  <- readRDS("test_paths_wk.rds")
PAR            <- readRDS("PAR.rds")
r              <- readRDS("r.rds")
RHO1           <- readRDS("RHO1.rds")
RHO2           <- readRDS("RHO2.rds")

par_est         <- PAR
par_est$mu1     <- PAR$mu1 + 0.05
par_est$lambda1 <- PAR$lambda1 * 1.35

# MISC 

wonham_step <- function(q1, q2, dBhat, lam1, lam2, r1, r2) {
  dr  <- r1 - r2
  q1n <- q1 + (lam2*exp(-q1) - (lam1+lam2) - 0.5*dr^2*(1-exp(q1))^2)*DT +
    dr*(1-exp(q1))*dBhat
  q2n <- q2 + ((lam1+lam2)*(exp(-q2)-1) - lam2*exp(-q2) -
                 0.5*dr^2*(1-exp(q2))^2)*DT -
    dr*(1-exp(q2))*dBhat
  p   <- exp(q1n)/(exp(q1n)+exp(q2n))
  list(q1 = q1n, q2 = q2n, p = p)
}

sim_path <- function(mu1, mu2, sigma, lam1, lam2) {
  reg          <- integer(K)
  cur          <- 1L
  switch_probs <- c(lam1, lam2) * DT
  for (i in seq_len(K)) {
    if (runif(1L) < switch_probs[cur]) cur <- 3L - cur
    reg[i] <- cur
  }
  mu_vec <- ifelse(reg == 1L, mu1, mu2)
  ret    <- (mu_vec - sigma^2/2)*DT + sigma*sqrt(DT)*rnorm(K)
  list(ret = ret, regime = reg)
}


lnf_th <- function(t, p, th) {
  tau <- T_HOR - t
  th[1]*tau + th[2]*tau^2 + th[3]*p*tau  + th[4]*p*tau^2 +
    th[5]*p^2*tau + th[6]*p^2*tau^2
}
f_th      <- function(t, p, th) exp(pmin(lnf_th(t, p, th), 500))
dlnf_dp   <- function(t, p, th) {          # fp/f = ∂lnf^θ/∂p
  tau <- T_HOR - t
  th[3]*tau + th[4]*tau^2 + 2*th[5]*p*tau + 2*th[6]*p*tau^2
}
# No change needed — both already handle vector t via element-wise ops.
# Explicitly confirm by adding a comment:
g_th <- function(t, th) { tau <- T_HOR-t; th[7]*tau + th[8]*tau^2 }  # vector-safe

V_th <- function(t, x, p, th, om)
  (x-om)^2 / f_th(t, p, th) + g_th(t, th) - (om-Z)^2               # vector-safe

# Gradient ∂V^θ/∂θ (equations 33-34)
dV_dth <- function(t, x, p, th, om) {
  tau <- T_HOR-t; ft <- f_th(t,p,th); cx <- -(x-om)^2/ft
  c(cx*c(tau, tau^2, p*tau, p*tau^2, p^2*tau, p^2*tau^2),   # θ_f part
    c(tau, tau^2))                                             # θ_g part
}

# Policy mean μ^φ, variance σ^φ
p_mean <- function(t, x, p, phi, th, om) {
  adj <- p - dlnf_dp(t,p,th)*p*(1-p)
  -(x-om)*(phi[1]*adj + phi[2])
}
p_var  <- function(t, p, phi, th) LAM*exp(phi[3])*f_th(t,p,th)
p_ent  <- function(t, p, phi, th) 0.5*log(2*pi*exp(1)*p_var(t,p,phi,th))
p_samp <- function(t, x, p, phi, th, om)
  rnorm(1L, p_mean(t,x,p,phi,th,om), sqrt(p_var(t,p,phi,th)))

oos_eval <- function(par_sim, par_filt, phi, th, om,
                     mode = c("poemv","emv"), n_test = 1000L) {
  mode <- match.arg(mode)
  
  if (any(!is.finite(phi)) || any(!is.finite(th)))
    stop("oos_eval: phi or th contains NaN/Inf — check training output")
  
  r1f  <- (par_filt$mu1-r)/par_filt$sigma
  r2f  <- (par_filt$mu2-r)/par_filt$sigma
  XT   <- numeric(n_test)
  
  for (trial in seq_len(n_test)) {
    sim <- sim_path(par_sim$mu1, par_sim$mu2, par_sim$sigma,
                    par_sim$lambda1, par_sim$lambda2)
    xk <- X0; pk <- 0.5; q1 <- q2 <- log(0.5)
    for (k in seq_len(K)) {
      tk <- (k-1L)*DT
      uk <- if (mode=="poemv") {
        p_samp(tk, xk, pk, phi, th, om)
      } else {
        sg <- pmax(LAM*exp(phi[3])*
                     exp(pmin(th[1]*(T_HOR-tk)+th[2]*(T_HOR-tk)^2, 500)),
                   1e-10)
        rnorm(1L, -(xk-om)*phi[1], sqrt(sg))
      }
      uk    <- pmax(pmin(uk, 1e4), -1e4)
      ret_k <- sim$ret[k]
      xk    <- xk + uk*(exp(ret_k)-1)
      if (!is.finite(xk)) xk <- X0
      if (mode=="poemv") {
        mhat <- par_filt$mu1*pk+par_filt$mu2*(1-pk)
        dBh  <- (ret_k-(mhat-par_filt$sigma^2/2)*DT)/par_filt$sigma
        dBh  <- pmax(pmin(dBh, 5*sqrt(DT)), -5*sqrt(DT))
        ws   <- wonham_step(q1,q2,dBh,par_filt$lambda1,par_filt$lambda2,r1f,r2f)
        q1<-ws$q1; q2<-ws$q2; pk<-ws$p
      }
    }
    XT[trial] <- xk
  }
  c(mean=mean(XT), var=var(XT), sharpe=(mean(XT)-X0)/sd(XT))
}

# =============================================================
# REAL-DATA EPISODE SAMPLER
# Replaces sim_path() — samples a random overlapping 10-year
# block from the real weekly return paths
# =============================================================

sample_real_episode <- function(paths_wk) {
  path_idx <- sample(length(paths_wk), 1L)
  ret_full <- as.numeric(paths_wk[[path_idx]])
  ret_full <- ret_full[is.finite(ret_full)]   # strip NAs — only change needed
  n <- length(ret_full)
  if (n < K) {
    ret_full <- sample(ret_full, K, replace=TRUE)
    return(list(ret=ret_full, regime=rep(1L, K)))
  }
  start <- sample(1L:(n-K+1L), 1L)
  ret   <- ret_full[start:(start+K-1L)]
  list(ret=ret, regime=rep(NA_integer_, K))
}

# =============================================================
# POEMV on real data — same algorithm, real episodes
# =============================================================

poemv_real <- function(par_filt, paths_wk,
                       n_iter    = 10000L,
                       n_batch   = 5L,
                       s_omega   = 50L,
                       eta_th    = 1e-12,
                       eta_phi   = 5e-5,
                       eta_omega = 5e-3,
                       verbose   = TRUE) {
  
  r1f <- (par_filt$mu1-r)/par_filt$sigma
  r2f <- (par_filt$mu2-r)/par_filt$sigma
  
  th  <- rep(0, 8L); phi <- rep(0, 3L); om <- Z; X_T_acc <- 0
  hist <- data.frame(phi1=NA_real_, phi2=NA_real_,
                     omega=NA_real_, X_T=NA_real_)[rep(1,n_iter),]
  
  for (iter in seq_len(n_iter)) {
    
    # ── Decaying learning rates ───────────────────────────
    decay     <- 1 / (1 + 1e-5 * iter)
    eta_th_t  <- eta_th  * decay
    eta_phi_t <- eta_phi * decay
    
    # ── Mini-batch accumulators ───────────────────────────
    dth_acc  <- rep(0, 8L)
    dphi_acc <- rep(0, 3L)
    XT_batch <- 0
    
    for (b in seq_len(n_batch)) {
      
      ep  <- sample_real_episode(paths_wk)
      ret <- ep$ret
      
      X  <- numeric(K+1L); X[1L] <- X0
      P  <- numeric(K+1L); P[1L] <- 0.5
      U_ <- H_ <- V_ <- numeric(K)
      q1 <- q2 <- log(0.5)
      
      for (k in seq_len(K)) {
        tk  <- (k-1L)*DT; xk <- X[k]; pk <- P[k]
        uk  <- p_samp(tk, xk, pk, phi, th, om)
        uk  <- pmax(pmin(uk, 1e4), -1e4)
        xn  <- xk + uk*(exp(ret[k])-1)
        if (!is.finite(xn)) xn <- xk
        X[k+1L] <- xn
        
        mhat <- par_filt$mu1*pk + par_filt$mu2*(1-pk)
        dBh  <- (ret[k] - (mhat - par_filt$sigma^2/2)*DT) / par_filt$sigma
        dBh  <- pmax(pmin(dBh, 5*sqrt(DT)), -5*sqrt(DT))
        ws   <- wonham_step(q1, q2, dBh, par_filt$lambda1,
                            par_filt$lambda2, r1f, r2f)
        q1 <- ws$q1; q2 <- ws$q2; P[k+1L] <- ws$p
        U_[k] <- uk
        H_[k] <- p_ent(tk, pk, phi, th)
        V_[k] <- V_th(tk, xk, pk, th, om)
      }
      
      XT <- X[K+1L]
      if (!is.finite(XT)) XT <- X0
      hT <- (XT-om)^2 - (om-Z)^2
      
      # PE gradient
      cumH    <- rev(cumsum(rev(H_))) * DT
      td_vec  <- hT - V_ - LAM*cumH
      tau_vec <- T_HOR - (seq_len(K)-1L)*DT
      px_vec  <- P[1L:K]
      cx_vec  <- -(X[1L:K]-om)^2 / pmax(f_th(tau_vec, px_vec, th), 1e-10)
      G_mat   <- cbind(
        cx_vec*tau_vec,          cx_vec*tau_vec^2,
        cx_vec*px_vec*tau_vec,   cx_vec*px_vec*tau_vec^2,
        cx_vec*px_vec^2*tau_vec, cx_vec*px_vec^2*tau_vec^2,
        tau_vec,                 tau_vec^2
      )
      dth_acc <- dth_acc + colSums(td_vec * G_mat) * DT
      
      # PG gradient
      tk_vec   <- (seq_len(K)-1L)*DT
      xk_vec   <- X[1L:K]; pk_vec <- P[1L:K]
      adj_vec  <- pk_vec - dlnf_dp(tk_vec, pk_vec, th)*pk_vec*(1-pk_vec)
      mu_p_vec <- -(xk_vec-om)*(phi[1]*adj_vec + phi[2])
      sg_p_vec <- pmax(LAM*exp(phi[3])*f_th(tk_vec, pk_vec, th), 1e-10)
      res_vec  <- U_ - mu_p_vec
      sc12_vec <- -(res_vec*(xk_vec-om))/sg_p_vec
      Vn_vec   <- c(V_th((seq_len(K-1L))*DT, X[2L:K], P[2L:K], th, om), hT)
      td2_vec  <- Vn_vec - V_ - LAM*H_*DT
      dphi_acc <- dphi_acc + c(
        sum(sc12_vec * adj_vec * td2_vec),
        sum(sc12_vec * td2_vec),
        sum((res_vec^2 - sg_p_vec)/(2*sg_p_vec)*td2_vec) - LAM*0.5*K*DT
      )
      
      XT_batch <- XT_batch + XT
      
    } # end batch loop
    
    # ── Apply averaged + clipped gradients ────────────────
    dth_acc  <- dth_acc  / n_batch
    dphi_acc <- dphi_acc / n_batch
    
    dth_acc  <- pmax(pmin(dth_acc,  1e2), -1e2)
    dphi_acc <- pmax(pmin(dphi_acc, 1e2), -1e2)
    
    th  <- th  + eta_th_t  * dth_acc
    phi <- phi - eta_phi_t * dphi_acc
    
    phi[1] <- pmax(pmin(phi[1], 20), -20)
    phi[2] <- pmax(pmin(phi[2], 20), -20)
    phi[3] <- pmax(pmin(phi[3], 10), -10)
    
    # ── Lagrange update ───────────────────────────────────
    X_T_acc <- X_T_acc + XT_batch / n_batch
    if (iter %% s_omega == 0L) {
      om <- om - eta_omega*(X_T_acc/s_omega - Z)
      X_T_acc <- 0
    }
    
    hist[iter,] <- list(phi[1], phi[2], om, XT_batch/n_batch)
    if (verbose && iter %% 2500L == 0L)
      cat(sprintf("  [%5d] φ₁=%7.3f  φ₂=%7.3f  ω=%.3f  XT=%.3f\n",
                  iter, phi[1], phi[2], om, XT_batch/n_batch))
  }
  list(th=th, phi=phi, omega=om, hist=hist)
}

# =============================================================
# EMV on real data
# =============================================================

emv_real <- function(paths_wk,
                     n_iter    = 10000L,
                     n_batch   = 5L,
                     s_omega   = 50L,
                     eta_th    = 1e-13,
                     eta_phi   = 5e-5,
                     eta_omega = 1e-4,
                     om_max    = 5.0,
                     verbose   = TRUE) {
  
  th <- rep(0,4); phi <- rep(0,3); om <- Z; X_T_acc <- 0
  hist <- data.frame(phi1=NA_real_, omega=NA_real_,
                     X_T=NA_real_)[rep(1,n_iter),]
  
  fe <- function(t,th) exp(pmin(th[1]*(T_HOR-t)+th[2]*(T_HOR-t)^2, 500))
  ge <- function(t,th) { th[3]*(T_HOR-t)+th[4]*(T_HOR-t)^2 }
  Ve <- function(t,x,th,om) {
    (x-om)^2/pmax(fe(t,th),1e-10) + ge(t,th) - (om-Z)^2
  }
  
  for (iter in seq_len(n_iter)) {
    
    # ── Decaying learning rates ───────────────────────────
    decay     <- 1 / (1 + 1e-5 * iter)
    eta_th_t  <- eta_th  * decay
    eta_phi_t <- eta_phi * decay
    
    # ── Mini-batch accumulators ───────────────────────────
    dth_acc  <- rep(0, 4L)
    dph_acc  <- rep(0, 3L)
    XT_batch <- 0
    
    for (b in seq_len(n_batch)) {
      
      ep  <- sample_real_episode(paths_wk)
      ret <- ep$ret
      X   <- numeric(K+1L); X[1L] <- X0
      Uv  <- Hv <- Vv <- numeric(K)
      
      for (k in seq_len(K)) {
        tk  <- (k-1L)*DT; xk <- X[k]
        sg  <- pmax(LAM*exp(phi[3])*fe(tk,th), 1e-10)
        uk  <- rnorm(1L, -(xk-om)*phi[1], sqrt(sg))
        uk  <- pmax(pmin(uk, 1e4), -1e4)
        xn  <- xk + uk*(exp(ret[k])-1)
        if (!is.finite(xn)) xn <- xk
        X[k+1L] <- xn
        Uv[k] <- uk
        Hv[k] <- 0.5*log(2*pi*exp(1)*sg)
        Vv[k] <- Ve(tk, xk, th, om)
      }
      
      XT <- X[K+1L]
      if (!is.finite(XT)) XT <- X0
      hT <- (XT-om)^2-(om-Z)^2
      
      # PE gradient
      cumH    <- rev(cumsum(rev(Hv)))*DT
      td_vec  <- hT - Vv - LAM*cumH
      tau_vec <- T_HOR-(seq_len(K)-1L)*DT
      ft_vec  <- pmax(fe(tau_vec, th), 1e-10)
      cx_vec  <- -(X[1L:K]-om)^2/ft_vec
      G_mat   <- cbind(cx_vec*tau_vec, cx_vec*tau_vec^2,
                       tau_vec,        tau_vec^2)
      dth_acc <- dth_acc + colSums(td_vec*G_mat)*DT
      
      # PG gradient
      sg_vec  <- pmax(LAM*exp(phi[3])*fe((seq_len(K)-1L)*DT,th), 1e-10)
      res_vec <- Uv-(-(X[1L:K]-om)*phi[1])
      sc1_vec <- -(res_vec*(X[1L:K]-om))/sg_vec
      Vn_vec  <- c(Ve((seq_len(K-1L))*DT, X[2L:K], th, om), hT)
      td2_vec <- Vn_vec-Vv-LAM*Hv*DT
      dph_acc <- dph_acc + c(
        sum(sc1_vec*td2_vec),
        0,
        sum((res_vec^2-sg_vec)/(2*sg_vec)*td2_vec) - LAM*0.5*K*DT
      )
      
      XT_batch <- XT_batch + XT
      
    } # end batch loop
    
    # ── Apply averaged + clipped gradients ────────────────
    dth_acc <- dth_acc / n_batch
    dph_acc <- dph_acc / n_batch
    
    dth_acc <- pmax(pmin(dth_acc, 1e2), -1e2)
    dph_acc <- pmax(pmin(dph_acc, 1e2), -1e2)
    
    th  <- th  + eta_th_t  * dth_acc
    phi <- phi - eta_phi_t * dph_acc
    
    phi[1] <- pmax(pmin(phi[1], 20), -20)
    phi[3] <- pmax(pmin(phi[3], 10), -10)
    
    # ── Lagrange update ───────────────────────────────────
    X_T_acc <- X_T_acc + XT_batch / n_batch
    if (iter %% s_omega == 0L) {
      om <- min(om - eta_omega*(X_T_acc/s_omega - Z), om_max)
      X_T_acc <- 0
    }
    
    hist[iter,] <- list(phi[1], om, XT_batch/n_batch)
    if (verbose && iter %% 2500L == 0L)
      cat(sprintf("  EMV [%5d] φ=%.3f  ω=%.3f  XT=%.3f\n",
                  iter, phi[1], om, XT_batch/n_batch))
  }
  list(th=th, phi=phi, omega=om, hist=hist)
}
## Call 
# =============================================================
# RUN — real data training
# =============================================================

cat("\n── Real Data: POEMV (true filter params) ───────\n")
res_T_real <- poemv_real(PAR, train_paths_wk, n_iter=100000L, n_batch=5L)
saveRDS(res_T_real, "res_T_real.rds")

cat("\n── Real Data: POEMV (estimated filter params) ──\n")
res_E_real <- poemv_real(par_est, train_paths_wk, n_iter=100000L, n_batch=5L)

saveRDS(res_E_real, "res_E_real.rds")

cat("\n── Real Data: EMV baseline ──────────────────────\n")
res_EMV_real <- emv_real(train_paths_wk, n_iter=100000L, n_batch=5L)

saveRDS(res_EMV_real, "res_EMV_real.rds")
# Results 

# =============================================================
# OOS EVALUATION — real-data trained policies on synthetic paths
# =============================================================
phi1_star <- (RHO1-RHO2)/PAR$sigma
phi2_star <- RHO2/PAR$sigma

cat("\nEvaluating out-of-sample on synthetic paths (real-data trained policies)…\n")
for (nm in c("res_T_real","res_E_real","res_EMV_real")) {
  obj <- get(nm)
  if (any(!is.finite(obj$phi)) || any(!is.finite(obj$th)))
    stop(sprintf("%s contains NaN/Inf — check training", nm))
}

oos_T_real   <- oos_eval(PAR,     PAR,
                         res_T_real$phi,   res_T_real$th,   res_T_real$omega,
                         "poemv", n_test=200000L)
oos_E_real   <- oos_eval(par_est, par_est,
                         res_E_real$phi,   res_E_real$th,   res_E_real$omega,
                         "poemv", n_test=200000L)
oos_EMV_real <- oos_eval(PAR,     PAR,
                         res_EMV_real$phi, res_EMV_real$th, res_EMV_real$omega,
                         "emv",   n_test=200000L)

cat("\n══ REAL DATA TABLE ═══════════════════════════════\n")
cat(sprintf("%-24s %7s %7s %7s %7s %7s %7s\n",
            "", "phi1","phi2","omega","Mean","Var","Sharpe"))
cat(sprintf("%-24s %7.3f %7.3f %7s %7s %7s %7s\n",
            "Theoretical target", phi1_star, phi2_star, "—","—","—","—"))
fmt <- "%-24s %7.3f %7.3f %7.3f %7.3f %7.3f %7.2f\n"
cat(sprintf(fmt,"POEMV (true filt)",
            res_T_real$phi[1],   res_T_real$phi[2],   res_T_real$omega,
            oos_T_real["mean"],   oos_T_real["var"],   oos_T_real["sharpe"]))
cat(sprintf(fmt,"POEMV (est filt)",
            res_E_real$phi[1],   res_E_real$phi[2],   res_E_real$omega,
            oos_E_real["mean"],   oos_E_real["var"],   oos_E_real["sharpe"]))
cat(sprintf("%-24s %7.3f %7s %7.3f %7.3f %7.3f %7.2f\n","EMV",
            res_EMV_real$phi[1], "—", res_EMV_real$omega,
            oos_EMV_real["mean"], oos_EMV_real["var"], oos_EMV_real["sharpe"]))
saveRDS(oos_T_real, "oos_T_real.rds")
saveRDS(oos_E_real, "oos_E_real.rds")
saveRDS(oos_EMV_real, "oos_EMV_real.rds")
# =============================================================
# CONVERGENCE PLOTS — real data
# =============================================================

roll_avg  <- function(x, w=500) stats::filter(x, rep(1/w,w), sides=2)
cols      <- c("steelblue","darkorange","darkgreen")
nms       <- c("POEMV (true)","POEMV (est)","EMV")
hst_list  <- list(res_T_real$hist, res_E_real$hist, res_EMV_real$hist)

par(mfrow=c(1,1), mar=c(4,4,3,1))

# φ₁ convergence
all_y <- unlist(lapply(hst_list, `[[`, "phi1"))
plot(NA, xlim=c(1, nrow(hst_list[[1]])),
     ylim=range(all_y, na.rm=TRUE),
     xlab="Iteration", ylab=expression(phi[1]),
     main=expression(phi[1] ~ "— Real Data Convergence"))
for (i in 1:3) {
  y <- hst_list[[i]]$phi1
  lines(y,           col=adjustcolor(cols[i], 0.15), lwd=0.5)
  lines(roll_avg(y), col=cols[i], lwd=2)
}
abline(h=phi1_star, lty=2, col="red")
legend("topleft",
       c("POEMV true (smoothed)","POEMV est (smoothed)",
         "EMV (smoothed)", expression(phi[1]^"*")),
       lty=c(1,1,1,2), col=c(cols,"red"), bty="n", cex=0.8)
# Terminal wealth
par(mfrow=c(1,1), mar=c(4,4,3,1))
wlth      <- lapply(hst_list, function(h) roll_avg(h$X_T))
ylim_wlth <- c(0, quantile(unlist(wlth), 0.995, na.rm=TRUE)*1.1)
plot(NA, xlim=c(1, length(wlth[[1]])), ylim=ylim_wlth,
     xlab="Iteration", ylab="Terminal wealth (rolling avg)",
     main="Terminal Wealth — Real Data, All Methods")
for (i in 1:3) lines(wlth[[i]], col=cols[i], lwd=2)
abline(h=Z, lty=2, col="red")
legend("bottomright",
       c(nms, paste0("Target Z=",Z)),
       lty=c(1,1,1,2), col=c(cols,"red"), bty="n", cex=0.8)
par(mfrow=c(1,1))

# In Time 
# =============================================================
# OOS PATH EVALUATOR — real test paths, full wealth trajectory
# =============================================================
oos_paths <- function(par_filt, phi, th, om,
                      mode = c("poemv", "emv"),
                      paths_wk) {
  mode <- match.arg(mode)
  if (any(!is.finite(phi)) || any(!is.finite(th)))
    stop("oos_paths: phi or th contains NaN/Inf — check training output")
  
  r1f     <- (par_filt$mu1-r)/par_filt$sigma
  r2f     <- (par_filt$mu2-r)/par_filt$sigma
  n_paths <- length(paths_wk)
  W       <- matrix(NA_real_, nrow=n_paths, ncol=K+1L)
  
  for (trial in seq_len(n_paths)) {
    ret_full <- as.numeric(paths_wk[[trial]])
    ret_full <- ret_full[is.finite(ret_full)]
    n        <- length(ret_full)
    
    if (n < K) {
      ret <- sample(ret_full, K, replace=TRUE)
    } else {
      ret <- ret_full[(n-K+1L):n]
    }
    
    xk <- X0; pk <- 0.5; q1 <- q2 <- log(0.5)
    W[trial, 1L] <- xk
    
    for (k in seq_len(K)) {
      tk <- (k-1L)*DT
      uk <- if (mode=="poemv") {
        p_samp(tk, xk, pk, phi, th, om)
      } else {
        sg <- pmax(LAM*exp(phi[3])*
                     exp(pmin(th[1]*(T_HOR-tk)+th[2]*(T_HOR-tk)^2, 500)),
                   1e-10)
        rnorm(1L, -(xk-om)*phi[1], sqrt(sg))
      }
      uk <- pmax(pmin(uk, 1e4), -1e4)
      xk <- xk + uk*(exp(ret[k])-1)
      if (!is.finite(xk)) xk <- W[trial, k]
      
      if (mode=="poemv") {
        mhat <- par_filt$mu1*pk + par_filt$mu2*(1-pk)
        dBh  <- (ret[k] - (mhat - par_filt$sigma^2/2)*DT) / par_filt$sigma
        dBh  <- pmax(pmin(dBh, 5*sqrt(DT)), -5*sqrt(DT))
        ws   <- wonham_step(q1, q2, dBh, par_filt$lambda1,
                            par_filt$lambda2, r1f, r2f)
        q1 <- ws$q1; q2 <- ws$q2; pk <- ws$p
      }
      W[trial, k+1L] <- xk
    }
  }
  W
}
# =============================================================
# RUN ON TEST PATHS
# =============================================================

cat("\nGenerating OOS wealth paths on real test data…\n")
paths_T   <- oos_paths(PAR, res_T_real$phi, res_T_real$th,
                       res_T_real$omega, "poemv", test_paths_wk)
paths_EMV <- oos_paths(PAR, res_EMV_real$phi, res_EMV_real$th,
                       res_EMV_real$omega, "emv", test_paths_wk)
saveRDS(paths_T, "paths_T.rds")
saveRDS(paths_EMV, "paths_EMV.rds")
# =============================================================
# FIGURE 7 ANALOG
# =============================================================

fig7 <- function(W_poemv, W_emv) {
  steps <- ncol(W_poemv)   # K+1
  
  mu_p  <- colMeans(W_poemv, na.rm=TRUE)
  sd_p  <- apply(W_poemv, 2, sd, na.rm=TRUE)
  mu_e  <- colMeans(W_emv,   na.rm=TRUE)
  sd_e  <- apply(W_emv,   2, sd, na.rm=TRUE)
  
  xs    <- seq_len(steps)
  ylim  <- range(c(mu_p - 0.5*sd_p, mu_p + 0.5*sd_p,
                   mu_e - 0.5*sd_e, mu_e + 0.5*sd_e), na.rm=TRUE)
  
  col_p <- "darkorange"; col_e <- "steelblue"
  
  plot(NA, xlim=c(1, steps), ylim=ylim,
       xlab="Weeks into out-of-sample horizon (2013–2022)",
       ylab="Portfolio wealth",
       main="Out-of-Sample Wealth: POEMV vs EMV")
  
  # EMV band
  polygon(c(xs, rev(xs)),
          c(mu_e + 0.5*sd_e, rev(mu_e - 0.5*sd_e)),
          col=adjustcolor(col_e, 0.2), border=NA)
  # POEMV band
  polygon(c(xs, rev(xs)),
          c(mu_p + 0.5*sd_p, rev(mu_p - 0.5*sd_p)),
          col=adjustcolor(col_p, 0.2), border=NA)
  
  lines(mu_e, col=col_e, lwd=2)
  lines(mu_p, col=col_p, lwd=2)
  abline(h=X0, lty=2, col="grey40")
  abline(h=Z,  lty=2, col="red")
  legend("topleft",
         c("POEMV (true)", "EMV", "Initial wealth", "Target z=2"),
         lty=c(1,1,2,2), lwd=c(2,2,1,1),
         col=c(col_p, col_e, "grey40", "red"), bty="n")
}

fig7(paths_T, paths_EMV)
fig7 <- function(W_poemv, W_emv) {
  steps <- ncol(W_poemv)   # K+1
  
  # Build weekly date axis from start of 2013
  date_start <- as.Date("2013-01-01")
  dates <- seq(date_start, by = "week", length.out = steps)
  xs <- as.numeric(dates)   # numeric for polygon/lines; axis handles labels
  
  mu_p  <- colMeans(W_poemv, na.rm = TRUE)
  sd_p  <- apply(W_poemv, 2, sd, na.rm = TRUE)
  mu_e  <- colMeans(W_emv,   na.rm = TRUE)
  sd_e  <- apply(W_emv,   2, sd, na.rm = TRUE)
  
  ylim <- range(c(mu_p - 0.5*sd_p, mu_p + 0.5*sd_p,
                  mu_e - 0.5*sd_e, mu_e + 0.5*sd_e), na.rm = TRUE)
  
  col_p <- "darkorange"; col_e <- "steelblue"
  
  plot(NA,
       xlim = range(xs), ylim = ylim,
       xlab = "Date", ylab = "Portfolio wealth",
       main = "Out-of-Sample Wealth: POEMV vs EMV (2013–2022)",
       xaxt = "n")   # suppress default numeric axis
  
  # Yearly tick marks
  year_breaks <- seq(as.Date("2013-01-01"), as.Date("2023-01-01"), by = "year")
  axis.Date(1,
            at     = year_breaks,
            format = "%Y",
            las    = 1)
  
  # EMV band
  polygon(c(xs, rev(xs)),
          c(mu_e + 0.5*sd_e, rev(mu_e - 0.5*sd_e)),
          col = adjustcolor(col_e, 0.2), border = NA)
  # POEMV band
  polygon(c(xs, rev(xs)),
          c(mu_p + 0.5*sd_p, rev(mu_p - 0.5*sd_p)),
          col = adjustcolor(col_p, 0.2), border = NA)
  
  lines(xs, mu_e, col = col_e, lwd = 2)
  lines(xs, mu_p, col = col_p, lwd = 2)
  
  abline(h = X0, lty = 2, col = "grey40")
  abline(h = Z,  lty = 2, col = "red")
  
  legend("topleft",
         c("POEMV (true)", "EMV", "Initial wealth", "Target z=2"),
         lty = c(1, 1, 2, 2), lwd = c(2, 2, 1, 1),
         col = c(col_p, col_e, "grey40", "red"), bty = "n")
}

fig7(paths_T, paths_EMV)