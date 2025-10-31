# getwd()
source("C:/Users/samue/OneDrive/Desktop/Iowa_State_PS/STAT 5200/Misc/basicglm.txt")
help(glm)
# set.seed(43)
# jxmat<-matrix(c(rep(1,50),seq(0.1,15,15/50)),50,2,byrow=F)
jxmat<-matrix(c(rep(1,35),seq(0.1,15,15/35)),35,2,byrow=F)
examp1<-simbasicglm(c(0.5,0.1),jxmat,2.0,8,5,pwr=0.25)
exampresult1<-basic.glm(jxmat,examp1,8,5,pwr=0.25)
eta<-0.5+0.1*jxmat[,2]
lines(jxmat[,2],true)
lines(jxmat[,2],exampresult1$vals$muhat,lty=2)
legend(1,15,legend=c("True","Estimated"),lty=c(1,2))

# jxmat[1:15,]
# examp1<-simbasicglm(c(0.5,0.1),jxmat,2.0,8,5,pwr=0.25)

# Poisson
# Log Link
# Try \beta_0 = 0.4
# Try \beta_1 = 0.025

j1 <- simbasicglm(c(0.4,0.05),jxmat,1,2,1)
# j1
plot(jxmat[,2], j1)

jlres <- basic.glm(jxmat, examp1, 1, 1)
jlres <- basic.glm(jxmat, examp1, 2, 1)

