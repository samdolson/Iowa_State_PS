source("C:/Users/samue/OneDrive/Desktop/Iowa_State_PS/STAT 5200/Misc/newtraph.txt")

dergamfctn<-function(ps){
  a<-ps[1]
  b<-ps[2]
  ga<-gamma(a)
  gb<-gamma(b)
  gab<-gamma(a+b)
  gprime<-function(x,par){ (x^{par-1})*(log(x))*exp(-x) }
  g2prime<-function(x,par){ (x^{par-1})*(log(x)^2)*exp(-x) }
  gpa<-integrate(gprime,0,Inf,par=a)$value
  gpb<-integrate(gprime,0,Inf,par=b)$value
  gpab<-integrate(gprime,0,Inf,par=(a+b))$value
  g2pa<-integrate(g2prime,0,Inf,par=a)$value
  g2pb<-integrate(g2prime,0,Inf,par=b)$value
  g2pab<-integrate(g2prime,0,Inf,par=(a+b))$value
  res1<-c(ga,gb,gab,gpa,gpb,gpab,g2pa,g2pb,g2pab)
  res<-matrix(res1,3,3,byrow=F)
  return(res)
}

betaders<-function(ps, y){
  n<-length(y)
  y1<-sum(log(y)); y2<-sum(log(1-y))
  iall<-dergamfctn(ps)
  gan<-iall[1,1]; gbn<-iall[2,1]; gabn<-iall[3,1]
  fa<-y1+n*((iall[3,2]/gabn)-(iall[1,2]/gan))
  fb<-y2+n*((iall[3,2]/gabn)-(iall[2,2]/gbn))
  fab<-n*(((iall[3,3]*gabn)-(iall[3,2]^2))/(gabn^2))
  faa<-fab-n*((((iall[1,3]*gan)-(iall[1,2]^2))/gan^2))
  fbb<-fab-n*((((iall[2,3]*gbn)-(iall[2,2]^2))/gbn^2))
  fm<-c(faa,fab,fab,fbb)
  dim(fm)<-c(2,2)
  logl<-((ps[1]-1)*y1)+((ps[2]-1)*y2)+n*(log(gabn)-log(gan)-log(gbn))
  res<-list(logl,c(fa,fb),fm)
  return(res)
}

set.seed(43) 

jdat<-rbeta(50,6,4)
mean(jdat)
var(jdat)

jout<-newtraph(betaders,jdat,c(1,1))

jout[[1]][1]-1.96*sqrt(jout[[3]][1,1])

jout[[1]][1]+1.96*sqrt(jout[[3]][1,1])

us<-seq(0.0001,0.9999,0.0001)
fs<-dbeta(us,jout[[1]][1],jout[[1]][2])
hist(jdat,prob=T,xlim=c(0,1),xlab="Variate Value",ylab="Density",main="")
lines(us,fs)

