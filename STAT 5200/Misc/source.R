weibullnegloglik<-function(pars,dat){
  #negative log likelihood for one sample Weibull model
  #for use with optim and nlm
  #pars is alpha, beta
  #dat is vector of responses
  #
  alp<-pars[1]; bet<-pars[2]
  lis<-alp*log(bet) + log(alp) + (alp-1)*log(dat) - (bet*dat)^alp
  llik<-sum(lis)
  nllik<-(-1)*llik
  return(nllik)
}
#----------------------------------------------------------------------
muint<-function(ahat,bhat,Invinf){
  #compute 95% interval estimate for Weibull expected value from one group
  #ahat, bhat are mles
  #Invinf is inverse information
  #
  aarg<-(ahat+1)/ahat
  mu<-(1/bhat)*gamma(aarg)
  cat("mu:",mu,fill=TRUE)
  dalp<-(1/bhat)*digamma(aarg)*gamma(aarg)*(-1/ahat^2)
  dbet<-(-1/bhat^2)*gamma(aarg)
  D<-c(dalp,dbet)
  vmu<-D%*%Invinf%*%D
  low<-mu-1.96*sqrt(vmu)
  up<-mu+1.96*sqrt(vmu)
  intvl<-c(low,up)
  return(intvl)
}
#-------------------------------------------------------------------------
survconfband<-function(ahat,bhat,Invinf,ts){
  #pointwise confidence band for Weibull survival function
  #ahat and bhat are estimates, ts are evaluation times
  #Infinf is inverse observed information
  #(not necessarily observed)
  #
  S<-exp(-(bhat*ts)^ahat)
  dalp<-S*(-(bhat*ts)^ahat)*log(bhat*ts)
  dbet<-S*(-ahat*(bhat*ts)^(ahat-1))*ts
  T<-length(ts)
  lows<-NULL; ups<-NULL
  cnt<-0
  repeat{
    cnt<-cnt+1
    tS<-S[cnt]
    tda<-dalp[cnt]
    tdb<-dbet[cnt]
    tD<-c(tda,tdb)
    tv<-tD%*%Invinf%*%tD
    tlow<-tS-1.96*sqrt(tv)
    tup<-tS+1.96*sqrt(tv)
    lows<-c(lows,tlow)
    ups<-c(ups,tup)
    if(cnt==T) break
  }
  res<-data.frame(S=S,low=lows,up=ups)
  return(res)
}
#-------------------------------------------------------------------------------
weibullders<-function(pars,dat){
  #compute log lik, first and second derivatives for use with
  #function newtraph
  #pars are alpha, beta
  #dat is vector of responses (times) for one group
  #
  alp<-pars[1]; bet<-pars[2]
  loglik<-(-1)*weibullnegloglik(pars,dat)
  dalps<-log(bet) + (1/alp) + log(dat) - ((bet*dat)^alp)*log(bet*dat)
  dbets<-(alp/bet)-alp*((bet*dat)^(alp-1))*dat
  d2alps<-(-1/alp^2) - ((bet*dat)^alp)*(log(bet*dat))^2
  d2bets<-(-alp/bet^2)-alp*(alp-1)*((bet*dat)^(alp-2))*dat^2
  d2alpbets<-(1/bet)-alp*((bet*dat)^(alp-1))*log(bet*dat)*dat-(1/bet)*(bet*dat)^alp
  dalp<-sum(dalps)
  dbet<-sum(dbets)
  d2alp<-sum(d2alps)
  d2bet<-sum(d2bets)
  d2alpbet<-sum(d2alpbets)
  grad<-c(dalp,dbet)
  H<-matrix(c(d2alp,d2alpbet,d2alpbet,d2bet),2,2,byrow=T)
  res<-list(loglik,grad,H)
  return(res)
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------


"newtraph" <- 
  function(ders, dat, x0)
  {
    cat("While N-R may be used for either minimization or\nmaximization")
    cat("The checks for progress in this function are written")
    cat("for maximization.  If you want to minimize, chage your")
    cat("derivative calculations (multiply by -1).")
    crit1<-1e-10
    crit2<-1e-06
    crit3<-1e-06
    c1<-0
    c2<-0
    c3<-0
    curnt <- x0
    nump<-length(x0)
    k <- 0
    repeat {
      k <- k + 1
      cat(" ",fill=T)
      cat(" ",fill=T)
      cat("Current estimates beginning iteration ", k, ":", fill = T)
      cat(curnt, fill = T)
      cat(" ",fill=T)
      int <- ders(curnt, dat)
      logL<-int[[1]]
      gi<-int[[2]]
      cat("Log likelihood for these estimates: ", fill = T)
      cat(logL, fill = T)
      cat("Gradient for these estimates: ",fill=T)
      cat(gi,fill=T)
      cat(" ", fill = T)
      Gi <- int[[3]]
      GiI <- solve(Gi)
      step <- GiI %*% gi
      new <- curnt - step
      sc <- 1
      repeat {
        sc <- sc + 1
        check <- ders(new, dat)
        if(check[[1]] <= logL) {new <- curnt - (1/sc) * step}
        if(check[[1]] > logL){newL<-check[[1]]
        newg<-check[[2]]}
        if(check[[1]] > logL) break
        if(sc == 10) {cat("Step halving not effective, try new starting values", fill = T)
          stop()}                        
      }
      dist1<-newL-logL
      dist2 <- (sum((new - curnt)^2))^0.5
      if(crit1 > dist1){c1<-1
      cat("Convergence criterion of ",crit1,"met for change in log likelihood",fill=T)}
      if(crit2 > dist2){c2<-1
      cat("Convergence criterion of ",crit2," met for change in estimates",fill=T)}
      if(sum(crit3 > abs(newg)) == nump){ c3<-1
      cat("Convergence criterion of ",crit3," met for sum of derivatives",fill=T)}
      if(c1+c2+c3==3)	break
      curnt <- new
    }
    cat("", fill = T)
    cat(" ",fill=T)
    final <- ders(new, dat)
    flogL<-final[[1]]
    fgrad<-final[[2]]
    fInf<--1*solve(final[[3]])
    cat("Final Estimates Are: ", new, fill = T)
    cat("", fill = T)
    cat("Final Log Likelihood: ", flogL, fill = T)
    cat("", fill = T)
    cat("Value of Gradient at Convergence:", fill = T)
    cat(fgrad, fill = T)
    cat("", fill = T)
    cat("Inverse Observed Information: ", fill = T)
    cat("(i.e., Inverse of Negative Hessian)", fill = T)
    cat("", fill = T)
    print(fInf)
    res<-list(new,flogL,fInf)
    return(res)
  }

