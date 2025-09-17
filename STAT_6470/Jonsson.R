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
