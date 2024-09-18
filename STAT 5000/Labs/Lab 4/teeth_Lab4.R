library(readr)
guinea_teeth <- read_csv("guinea_teeth.csv", col_types = cols(trt = 
                                   col_factor(levels = c("0", "1"))))
View(guinea_teeth)

sd(guinea_teeth$growth[guinea_teeth$trt=="1"])/
  sd(guinea_teeth$growth[guinea_teeth$trt=="0"])

var.test(x=guinea_teeth$growth[guinea_teeth$trt=="1"],
         y=guinea_teeth$growth[guinea_teeth$trt=="0"], 
         alternative="greater")

BF.var.test <- function(dat.response, dat.treatment){
  n1 = length(dat.response[dat.treatment==levels(dat.treatment)[1]])
  n2 = length(dat.response[dat.treatment==levels(dat.treatment)[2]])
  M = c(rep(median(dat.response[dat.treatment==levels(dat.treatment)[1]]), n1),
        rep(median(dat.response[dat.treatment==levels(dat.treatment)[2]]),n2))
  Z = abs(c(dat.response[dat.treatment==levels(dat.treatment)[1]], 
            dat.response[dat.treatment==levels(dat.treatment)[2]]) - M)
  G = c(dat.treatment[dat.treatment==levels(dat.treatment)[1]], 
        dat.treatment[dat.treatment==levels(dat.treatment)[2]])
  df = length(Z)-2
  BFstat = (t.test(Z~G, var.equal=T)$statistic)^2
  pval = pf(BFstat, 1, df, lower.tail=F)
  return(data.frame(BFstat=BFstat, pval=pval, row.names="results:"))
}
BF.var.test(guinea_teeth$growth, guinea_teeth$trt)

boxplot(guinea_teeth$growth ~ guinea_teeth$trt, xlab="Treatment", 
        ylab="Growth", main="Guinea Pig Teeth Experiment")
par(mfrow=c(2,1))
hist(guinea_teeth$growth[guinea_teeth$trt=="1"],
     main="Treatment = orange juice", xlab="Growth")
hist(guinea_teeth$growth[guinea_teeth$trt=="0"],
     main="Treatment = ascorbic acid", xlab="Growth")
par(mfrow=c(1,2))
qqnorm(scale(guinea_teeth$growth[guinea_teeth$trt=="1"]),
       main="Treatment = orange juice")
abline(a=0, b=1, col="red")
qqnorm(scale(guinea_teeth$growth[guinea_teeth$trt=="0"]),
       main="Treatment = ascorbic acid")
abline(a=0, b=1, col="red")

library(tidyverse)
library(moments)
numerical_stats <- guinea_teeth |> 
  group_by(trt) |>
  summarize(
    Y_mean = mean(growth),
    Y_med = quantile(growth, 0.5),
    Y_sd = sd(growth),
    Y_IQR = quantile(growth, 0.75) - quantile(growth, 0.25),
    Y_skew = skewness(growth),
    Y_kurt = kurtosis(growth),
    Y_excess = Y_kurt - 3
  )
numerical_stats

shapiro.test(guinea_teeth$growth[guinea_teeth$trt=="1"])
shapiro.test(guinea_teeth$growth[guinea_teeth$trt=="0"])


t.test(growth~trt, data=guinea_teeth, var.equal=T)

t.test(growth~trt, data=guinea_teeth, var.equal=F)

wilcox.test(growth~trt, data=guinea_teeth, exact=F)
# sum of the ranks in "first" group is actually W + n*(n+1)/2
19.5+(10*11/2)

