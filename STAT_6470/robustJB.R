ft.RJB.multi <- function(x)
{
    xm <- scale(x, scale = FALSE)
    x.pc <- prcomp(xm, retx = T)
    values <- x.pc$sd
    z <-  x.pc$x %*% diag(1/sqrt(values)) 
    result <- apply(z, 1, ft.robustJB)
    stat <- sum(result)
    stat 
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
