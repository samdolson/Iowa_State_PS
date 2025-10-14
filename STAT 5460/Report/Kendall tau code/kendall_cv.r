#An adaptive test based on Kendall's tau for independence in high dimensions
rm(list=ls())
set.seed(123456)
library(RcppArmadillo)
library(Rcpp)
library(MASS)
setwd("D:/R/kenall tau") 
Rcpp::sourceCpp("taucpp.cpp")
t1 <- Sys.time()
RES <- NULL
N <- c(50,100)
P <- c(50,100,200,400)
M <- 5000
alpha <- 0.05
for(n in N){
  REs <- NULL
  for(p in P){ 
    
    ECP <- matrix(0,M,3)
    Pnp<-rep(0,2)
    Tn <- rep(0,2)
    Fn <- matrix(0,M,3)
    for(m in 1:M){
      cat(n,p,m,"\r")
      X <- matrix(rnorm(n*p),n,p) 
      Tn <- taucpp(X)
      Fn[m,1]<-Tn[1]
      Fn[m,2] <- Tn[2]-4*log(p)+log(log(p))
      Pnp[1] <- 1-pnorm(Tn[1])
      Pnp[2] <- 1-exp(-exp(-(Tn[2]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
      Cnp=min(Pnp[1],Pnp[2])
      Fn[m,3]<-Cnp
      

      qn <- c(quantile(Fn[,1], c(0.025,0.975)),quantile(Fn[,2], 0.95),quantile(Fn[,3],  0.05))
      
    }
    REs <- rbind(REs,qn)
  }
  rownames(REs) <- paste("p=",P)
  RES <- rbind(RES,REs)
}

t2 <- Sys.time()
print((t2-t1))
print(RES)

write.csv(RES,file = "D:/R/kenall tau/data/cv.csv")





