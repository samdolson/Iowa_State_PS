library(readr) 
library(dplyr)
source("C:/Users/samue/OneDrive/Desktop/Iowa_State_PS/STAT 5200/Misc/source.R")

dataWei <- read.table("C:/Users/samue/OneDrive/Desktop/Iowa_State_PS/STAT 5200/Misc/weibulldat.txt", header = TRUE)

dat1 <- dataWei |> 
  filter(trt == 1)

dat2 <- dataWei |>
  filter(trt == 2)

stem(dat1$time)
stem(dat2$time)

result1 <- nlm(f = weibullnegloglik, p = c(1, 0.6), hessian = T, dat = dataWei)

result2 <- nlm(f = weibullnegloglik, p = c(1, 0.6), hessian = T, dat = dat1)

result3 <- nlm(f = weibullnegloglik, p = c(1, 0.6), hessian = T, dat = dat2)

# inverse observed information of negative log likelihood
# we have the negative log likelihood, just need to invert 

result2$estimate
result3$estimate 

InvertInf1 <- solve(result2$hessian)
InvertInf1

InvertInf2 <- solve(result3$hessian)
InvertInf2

# Validate Some of our Results 

result3 <- nlm(f = weibullnegloglik, p = c(1, 0.6), hessian = T, dat = dat2)
result3_1 <- nlm(f = weibullnegloglik, p = c(2, 1), hessian = T, dat = dat2)
result3_2 <- nlm(f = weibullnegloglik, p = c(0.1, 0.1), hessian = T, dat = dat2)

# check whether the iteratively derived mle is consistent (we hope it is!)
result3$code
result3$estimate

result3_1$code
result3_1$estimate

result3_2$code
result3_2$estimate

# Let's get a Wald-based interval 

result2$estimate - 1.96 * sqrt(diag(InvertInf1))
result2$estimate + 1.96 * sqrt(diag(InvertInf1))

result3$estimate - 1.96 * sqrt(diag(InvertInf2))
result3$estimate + 1.96 * sqrt(diag(InvertInf2))

# Now, want Likelihood Ratio Test B/W the two groups 
# I don't know what was happening here; I need to revisit when the pdf is available for comparison 

result2$minimum

result3$minimum

llikeF <- - result1$minimum - result2$minimum

llikeF

# Note: 
# can re-run the above using the other "black box" optimizer, `optim`; again, hope/expect the results to be the same 
optimRes1 <- optim(par = c(1, 0.6), fn = weibullnegloglik, hessian = T, dat = dat1)

optimRes1$par
result2$estimate

optimRes1$hessian
result2$hessian