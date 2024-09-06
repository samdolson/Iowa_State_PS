# observations.pdf
# used for todays lecture 

# for the data.frame example
# did not have constant variance of the error term
# not identically distributed (iid) 

rm()

x <- 1:20 
w <- 1 + sqrt(x)/2

dummy <- data.frame(x = x, y = x + rnorm(x)*w) 

fm <- lm(formula = y ~ x, data = dummy) 
summary(fm) 

plot(dummy)

residuals(fm)
plot(residuals(fm))

attach(dummy)

plot(x = x, y = y)
abline(a = 0, b = 1, lty = 3)

plot(x = x, y = w)

detach()

plot(x = fitted(fm), y = resid(fm), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted")
# abline(a = 0, b = 0)
abline(h = 0)

rm()

# rm(fm, x, y, dummy)

# Review the idea of a "Structure" in C (C++) 
changeRes <- rep(NA, length(res))

for(i in 1:length(changeRes)) {
  changeRes[i] <- res[i] - res[i+1]
}

newChange <- na.omit(changeRes)


# reading data
setwd("C:/Users/samue/OneDrive/Desktop/Iowa_State_PS/STAT 5970/Notes")

Houses <-read.table(file="houses.dat", header = T)
# HousesOnline <- read.table(file = "WEBADDRESS", header = T)

