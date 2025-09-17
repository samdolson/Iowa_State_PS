# 647 Normality Tests R script (one call)

#
# we will now do an example applying Cramer's characterization of the normal
# distribution, but we will do this in some sort of efficient manner so that R
# will be happy.
#
# we generate N = 100,000 random unit projections in p-D
#

testnormality <- function(X, numproj = 5e2, parallel = TRUE, ncores = NULL) {
  # note that the value returned is the q-value of the test
  # alternative to Shapiro-Wilks, we could also use energy test for 1D, but that is not implemented. 
  p <- ncol(X)
  n <- nrow(X)
  
  if (parallel) {    
    require(parallel)
    if (is.null(ncores))
      ncores <- detectCores()
    clus <- makeCluster(ncores)
    x <- parApply(cl = clus, X = as.matrix(rep(numproj, p)), MARGIN = 1, FUN = function(x) (rnorm(n = x)))
    # generate numproj standard
    # p-variate
    # normal random variables.
    #        dim(x)  <- rev(dim(x)) #change so that the dimensions are what we want for the next step.
    y <- sqrt(parApply(cl = clus, X = x^2, MARGIN = 1, FUN = sum))
    # not completely sure if rowSums is not faster
    z <- x / y
    tempdat <- z %*% t(as.matrix(X))  ## this gives rise to a numproj x p
    ## matrix called tempdat here. we now
    ## perform Shapiro-Wilks' test and calculate individual p-values on
    ## each of numproj observation sets.
    pvals <- parApply(cl = clus, X = tempdat, MARGIN = 1, FUN = function(x)(shapiro.test(x)$p.value))
  } else {
    x <- matrix(rnorm(numproj * p), ncol = p)
    y <- sqrt(rowSums(x^2))
    z <- x / y
    tempdat <- z %*% t(as.matrix(X)) 
    pvals <- apply(X = tempdat, MARGIN = 1, FUN = function(x) (shapiro.test(x)$p.value))
  }
  min(p.adjust(pvals, method = "holm"))
}



##
## Generalized Shapiro-Wilk's test.
##
tr.sw <- function(x) {
  sw <- shapiro.test(x = x)
  c(tstat = sw$statistic, p.value = sw$p.value)
}

gen.SW <-  function(X)
{
  library(Compositional)
  n <- dim(X)[1]
  p <- dim(X)[2]
  A <- helm(n)
  Y <- A%*%X 
  Y2 <- crossprod(Y)/(n-1)
  ei <- eigen(Y2, symmetric = T)
  h <- ei$values #already in desc order
  d <- ei$vectors[, which(h > 0)] #Vectors are in the columns
  
  Zi <- Y%*%d
  
  SW.pval <- as.data.frame(t(apply(X = Zi, MARGIN = 2, FUN = tr.sw)))
  SW.pvalmin <- which.min(SW.pval$p.value)
  list(GSW.stat = SW.pval[SW.pvalmin,1], p.value =  SW.pval[SW.pvalmin,2])
}



p.val.JB <- function(x) {
  if ((NCOL(x) > 1) || is.data.frame(x)) 
    stop("x is not a vector or univariate time series")
  if (any(is.na(x))) 
    stop("NAs in x")
  DNAME <- deparse(substitute(x))
  n <- length(x)
  m1 <- sum(x)/n
  m2 <- sum((x - m1)^2)/n
  m3 <- sum((x - m1)^3)/n
  m4 <- sum((x - m1)^4)/n
  b1 <- (m3/m2^(3/2))^2
  b2 <- (m4/m2^2)
  STATISTIC <- n * b1/6 + n * (b2 - 3)^2/24
  PVAL <- 1 - pchisq(STATISTIC, df = 2)
  PARAMETER <- 2
  METHOD <- "Jarque Bera Test"
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  PVAL
}


Jonsson <- function(x) {
  JB.p.values <- apply(x, 2, p.val.JB)
  LM.stat <- -2*sum(log(JB.p.values))  
  pv <- 1 - pchisq(LM.stat, df = 2*dim(x)[2])
  pv
}  


ft.RJB.multi <- function(x, p)
{
  xm <- scale(x, scale = FALSE)
  x.pc <- prcomp(xm, retx = T)
  values <- x.pc$sd
  z <-  x.pc$x %*% diag(1/values) # changed based on class discussion 
  result <- apply(z, 2, ft.robustJB)
  stat <- sum(result)
  stat
  p.value <- 1-pchisq(stat, df = 2*p)
  p.value
}

ft.MAAD <- function(x)
{
  cc <- sqrt(pi/2)
  jj <- cc/length(x) * sum(abs(x - median(x)))
  jj
}

ft.robustJB <- function(x)
{
  nn <- length(x)
  mu3 <- sum((x - mean(x))^3)/nn
  mu4 <- sum((x - mean(x))^4)/nn
  jj <- ft.MAAD(x)
  RJB <- nn/6 * (mu3/jj^3)^2 + nn/64 * (mu4/jj^4 - 3)^2
  RJB
}



###############################################################
# Extreme and non-extreme BHEP multivariate test statistics
# h = smoothing parameter
# Mahalanobis = Mahalanobis distances and angles = D_ij
# n = sample size
# d = data-dimension

BHEP.infty = function(Mahalanobis,n,d)
{ 
  M.distances <- diag(Mahalanobis)
  b1 <- sum(Mahalanobis^3)/n^2
  b1t <- t(M.distances)%*%Mahalanobis%*%M.distances/n^2
  Tn0 <- b1/6 + b1t/4
  ##as.real(n*Tn0) #changed by RM
  as.double(n*Tn0)
}

BHEP.zero = function(Mahalanobis,n,d)
{
  Tninfty <- mean(exp(-diag(Mahalanobis)/2))
  sqrt(n)*abs(Tninfty - 2^(-d/2))
}

BHEP.stat = function(h,Mahalanobis,n,d) #corrected value
{
  M.distances <- diag(Mahalanobis)
  A <- matrix(rep(M.distances,n),ncol=n)
  M.differences <- A + t(A) - 2*Mahalanobis  #=D_ii-2D_ij+D_jj=||Y_i - Y_j||^2
  h1 <- 2*h^2
  h2 <- 2*h^2+1
  h3 <- 2*h^2+2
  term1 <- 2*sum(exp(-M.differences[upper.tri(M.differences)]/(2*h1)))/(n*h1^(d/2))
  term2 <- 2*sum(exp(-M.distances/(2*h2)))/(h2^(d/2))
  bhep <- (2*pi)^(d/2)*( term1 - term2 + h1^(-d/2) + n*h3^(-d/2) )
  bhep
}


###################################################################
# Multiple test P-value evaluation
# y = observed sample

p.value.BB = function(y)
{
  #sample size
  n <- length(y[,1])
  if (n<20) warning("sample size must be >= 20") 
  else {
    
    #data dimension
    d <- length(y[1,])
    if ((d==1)||(d>10)) warning("data dimension is not between 2 and 10") 
    else {
      
      #considered set of sizes
      alphas <- c(seq(0.0005,0.15,by=0.0005),seq(0.15+0.05,0.95,by=0.05))
      lalphas <- length(alphas)
      
      #considered set of sample sizes 
      ns <- c(seq(20,100,by=5),seq(110,200,by=10),seq(220,300,by=20),seq(350,800,by=50),900,1000)
      lns <- length(ns)
      
      line <- length(ns[ns<=n])
      n.close <- ns[line]
      
      if ((n.close==n)||(n>1000)) 
      {
        corr1 <- as.matrix(read.table("BZ_BB.R",skip=(d-2)*lns+line-1,nrows=1))
        corr2 <- as.matrix(read.table("BI_BB.R",skip=(d-2)*lns+line-1,nrows=1))
        corr3 <- as.matrix(read.table("BS_BB.R",skip=(d-2)*lns+line-1,nrows=1))
        corr4 <- as.matrix(read.table("BL_BB.R",skip=(d-2)*lns+line-1,nrows=1))
      }
      
      if ((n.close<n)&(n<=1000)) 
      {
        distanc<-(n-ns[line])/(ns[line+1]-ns[line]) #weighted distance between n and n.close
        corr1o <- as.matrix(read.table("BZ_BB.R",skip=(d-2)*lns+line-1,nrows=2)) 
        corr2o <- as.matrix(read.table("BI_BB.R",skip=(d-2)*lns+line-1,nrows=2))
        corr3o <- as.matrix(read.table("BS_BB.R",skip=(d-2)*lns+line-1,nrows=2))
        corr4o <- as.matrix(read.table("BL_BB.R",skip=(d-2)*lns+line-1,nrows=2))
        corr1 <- (1-distanc)*corr1o[1,]+distanc*corr1o[2,]
        corr2 <- (1-distanc)*corr2o[1,]+distanc*corr2o[2,]
        corr3 <- (1-distanc)*corr3o[1,]+distanc*corr3o[2,]
        corr4 <- (1-distanc)*corr4o[1,]+distanc*corr4o[2,]
      }
      
      #Test statistic evaluation
      Tn <- array(dim=c(lalphas))
      invS <- solve((1-1/n)*cov(y))
      med <- matrix(rep(colMeans(y),n),ncol=d,byrow=TRUE)
      Mahalanobis <- (y-med)%*%invS%*%t(y-med)
      s2 <- 3^(-d/2)-2^(-d)-d*2^(-d-3)
      
      Tn.a<-array(dim=c(4,lalphas))
      Tn.a[1,] <- BHEP.zero(Mahalanobis,n,d)^2/s2 - corr1
      Tn.a[2,] <- BHEP.infty(Mahalanobis,n,d) - corr2
      
      Tn.a[3,] <- BHEP.stat(0.448 + 0.026*d,Mahalanobis,n,d) - corr3
      Tn.a[4,] <- BHEP.stat(0.928 + 0.049*d,Mahalanobis,n,d) - corr4
      Tn <- apply(Tn.a,MARGIN=2,FUN=max)
      
      comp <- length(Tn[Tn<=0])
      
     # if ((comp==0)|(comp==1)) pv <- paste("p-value <",alphas[1]) else
     #   if (comp==lalphas) pv <- paste("p-value >=",alphas[comp]) else 
      #    pv <- paste(alphas[comp],"< p-value <",alphas[comp+1]) <- dumb
      
      if ((comp==0)|(comp==1)) pv <- alphas[1] else
        if (comp==lalphas) pv <- alphas[comp] else 
          pv <- alphas[comp+1]
      
      pv
      
    }}
}


###############################################################
# Extreme and non-extreme BHEP multivariate test statistics
# h = smoothing parameter
# Mahalanobis = Mahalanobis distances and angles = D_ij
# n = sample size
# d = data-dimension

BHEP.stat = function(h,Mahalanobis,n,d) #corrected value
{
  M.distances <- diag(Mahalanobis)
  A <- matrix(rep(M.distances,n),ncol=n)
  M.differences <- A + t(A) - 2*Mahalanobis  #=D_ii-2D_ij+D_jj=||Y_i - Y_j||^2
  h1 <- 2*h^2
  h2 <- 2*h^2+1
  h3 <- 2*h^2+2
  term1 <- 2*sum(exp(-M.differences[upper.tri(M.differences)]/(2*h1)))/(n*h1^(d/2))
  term2 <- 2*sum(exp(-M.distances/(2*h2)))/(h2^(d/2))
  bhep <- (2*pi)^(d/2)*( term1 - term2 + h1^(-d/2) + n*h3^(-d/2) )
  bhep
}


###############################################################
# Mardia's multivariate test statistics 
# Mahalanobis = Mahalanobis distances and angles = D_ij
# n = sample size
# d = data-dimension

MS.stat = function(Mahalanobis,n,d)
{
  sum(Mahalanobis^3)/n
}

MK.stat = function(Mahalanobis,n,d)
{
  sqrt(n)*abs(sum(diag(Mahalanobis)^2)/n - d*(d+2))
}


###################################################################
# Multiple test P-value evaluation
# y = observed sample

p.value.MB = function(y)
{
  #sample size
  n <- length(y[,1])
  if (n<20) warning("sample size must be >= 20") 
  else {
    
    #data dimension
    d <- length(y[1,])
    if ((d==1)||(d>10)) warning("data dimension is not between 2 and 10") 
    else {
      
      #considered set of sizes
      alphas <- c(seq(0.0005,0.15,by=0.0005),seq(0.15+0.05,0.95,by=0.05))
      lalphas <- length(alphas)
      
      #considered set of sample sizes 
      ns <- c(seq(20,100,by=5),seq(110,200,by=10),seq(220,300,by=20),seq(350,800,by=50),900,1000)
      lns <- length(ns)
      
      line <- length(ns[ns<=n])
      n.close <- ns[line]
      
      if ((n.close==n)||(n>1000)) 
      {
        corr1 <- as.matrix(read.table("MS_MB.R",skip=(d-2)*lns+line-1,nrows=1))
        corr2 <- as.matrix(read.table("MK_MB.R",skip=(d-2)*lns+line-1,nrows=1))
        corr3 <- as.matrix(read.table("BS_MB.R",skip=(d-2)*lns+line-1,nrows=1))
        corr4 <- as.matrix(read.table("BL_MB.R",skip=(d-2)*lns+line-1,nrows=1))
      }
      
      if ((n.close<n)&(n<=1000)) 
      {
        distanc<-(n-ns[line])/(ns[line+1]-ns[line]) #weighted distance between n and n.close
        corr1o <- as.matrix(read.table("MS_MB.R",skip=(d-2)*lns+line-1,nrows=2)) 
        corr2o <- as.matrix(read.table("MK_MB.R",skip=(d-2)*lns+line-1,nrows=2))
        corr3o <- as.matrix(read.table("BS_MB.R",skip=(d-2)*lns+line-1,nrows=2))
        corr4o <- as.matrix(read.table("BL_MB.R",skip=(d-2)*lns+line-1,nrows=2))
        corr1 <- (1-distanc)*corr1o[1,]+distanc*corr1o[2,]
        corr2 <- (1-distanc)*corr2o[1,]+distanc*corr2o[2,]
        corr3 <- (1-distanc)*corr3o[1,]+distanc*corr3o[2,]
        corr4 <- (1-distanc)*corr4o[1,]+distanc*corr4o[2,]
      }
      
      #Test statistic evaluation
      Tn <- array(dim=c(lalphas))
      invS <- solve((1-1/n)*cov(y))
      med <- matrix(rep(colMeans(y),n),ncol=d,byrow=TRUE)
      Mahalanobis <- (y-med)%*%invS%*%t(y-med)
      
      Tn.a<-array(dim=c(4,lalphas))
      Tn.a[1,] <- MS.stat(Mahalanobis,n,d) - corr1
      Tn.a[2,] <- MK.stat(Mahalanobis,n,d)^2 - corr2
      Tn.a[3,] <- BHEP.stat(0.448 + 0.026*d,Mahalanobis,n,d) - corr3
      Tn.a[4,] <- BHEP.stat(0.928 + 0.049*d,Mahalanobis,n,d) - corr4
      Tn <- apply(Tn.a,MARGIN=2,FUN=max)
      
      comp <- length(Tn[Tn<=0])
      
      if ((comp==0)|(comp==1)) pv <- alphas[1] else
        if (comp==lalphas) pv <- alphas[comp] else 
          pv <- alphas[comp+1]
      
      pv
      
    }}
}



