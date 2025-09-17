#
# we will now do an example applying Cramer's characterization of the normal
# distribution, but we will do this in some sort of efficient manner so that R
# will be happy.
#
# we generate N = 100,000 random unit projections in p-D
#

testnormality <- function(X, numproj = 1e5, parallel = TRUE, ncores = NULL) {
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
