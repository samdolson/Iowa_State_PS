# R Start Guide Example

library(readr)
salaries <- read_csv("C:/Users/samue/Downloads/salaries.csv")
View(salaries)

class(salaries)

# Change class of Sex from numeric to factor
class(salaries$Sex) = "Factor"

salaries$Salary
salaries$Sex

head(salaries)

mean(salaries$Salary[salaries$Sex==1])
mean(salaries$Salary[salaries$Sex==0])



hist(salaries$Salary)

par(mfrow = c(2,1)) 

hist(salaries$Salary[salaries$Sex==1], 
     xlab="dollars", main="Male's Salaries")

hist(salaries$Salary[salaries$Sex == 0],
     xlab="dollars", main="Female's Salaries")
