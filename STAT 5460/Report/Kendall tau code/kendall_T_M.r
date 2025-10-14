# An adaptive test based on Kendall's tau for independence in high dimensions
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
cv1<-c(-1.861, -1.854, -1.910, -1.885, -1.858, -1.948, -1.916, -1.913)
cv2<-c(2.006, 2.001, 1.971, 2.008, 1.975, 1.944, 1.945, 1.997)
cv3<-c(1.770, 1.242, 0.953, 0.798, 2.024, 1.963, 1.782, 1.679)
cv4<-c(0.029, 0.035, 0.038, 0.035, 0.031, 0.032, 0.037, 0.031)
Sim <- 2000
alpha <- 0.05
m=1
for(n in N){
  REs <- NULL
  for(p in P){ 
    
    ECP <- matrix(0,Sim,6)
    Pnp<-rep(0,2)
    Tn <- rep(0,2)
    for(sim in 1:Sim){
      cat(n,p,sim,"\r")
      #norm
     # X <- matrix(rnorm(n*p),n,p)
      
      #cauchy
      X<-matrix(rcauchy(n*p,0,1),n,p)
      
      #t(4)
      # X<-matrix(rt(n*p,4),n,p)
      
      #dnorm
      # Sigma <- matrix(0,p,p)
      # 
      # Sigma<-(1-log(p)/n/2)*diag(p)+log(p)/n/2*matrix(1,p,p)
      # 
      # X <- mvrnorm(n, rep(0,p),Sigma)
      
      
      #dcauchy
      # Z<-matrix(rcauchy(n*(p+1),0,1),n,(p+1))#cauchy(0,1)
      # X<-Z[,(1:p)]+2*log(p)/n*Z[,(2:(p+1))]
      
      #snorm
      # Sigma<-matrix(0,p,p)
      # Sigma<-1*diag(p)
      # Sigma[1,2]<-log(p)/8
      # Sigma[2,1]<-log(p)/8
      # X <- mvrnorm(n, rep(0,p),Sigma)
      
      
      #scauchy
      # Z<-matrix(rcauchy(n*p,0,1),n,p)#cauchy(0,1)
      # X<-cbind(Z[,1]+3*log(p)/n*Z[,2],Z[,2]+3*log(p)/n*Z[,1],Z[,3:p])
      
      Tn <- taucpp(X)
      ECP[sim,1] <- 1*(abs(Tn[1])>qnorm(1-alpha/2))
      ECP[sim,2] <- 1*((Tn[1]<cv1[m])||(Tn[1]>cv2[m]))
      
      ECP[sim,3] <- 1*((Tn[2]-4*log(p)+log(log(p)))>(-log(8*pi)-2*log(log(1/(1-alpha)))))
      ECP[sim,4] <- 1*((Tn[2]-4*log(p)+log(log(p)))>cv3[m])
      
      Pnp[1] <- 1-pnorm(Tn[1])
      Pnp[2] <- 1-exp(-exp(-(Tn[2]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
      Cnp=min(Pnp[1],Pnp[2])
      ECP[sim,5]<-1*(Cnp<1-sqrt(1-alpha))
      ECP[sim,6]<-1*(Cnp<cv4[m])
      
      ECP[sim,]<-c(ECP[sim,1],ECP[sim,2],ECP[sim,3],ECP[sim,4],ECP[sim,5],ECP[sim,6])
    }
    REs <- rbind(REs,colMeans(ECP))
    m=m+1;
  }
  rownames(REs) <- paste("p=",P)
  RES <- rbind(RES,REs)
}
colnames(RES) <- c("TS","MS","TM","MM","TC","MC")

t2 <- Sys.time()
print((t2-t1))
print(RES)

write.csv(t(RES),file = "D:/R/kenall tau/data/cauchy.csv")


