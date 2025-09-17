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
