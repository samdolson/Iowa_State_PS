library(readr)
hypertension <- read_csv("STAT 5000/Labs/Lab 1/Optional/hypertension.csv", 
                         col_types = cols(sodiumdiet = col_factor(levels = c("low", 
                                                                             "high"))))
View(hypertension)
print(hypertension)

summary(hypertension$bloodpressure[hypertension$sodiumdiet=="low"])
sd(hypertension$bloodpressure[hypertension$sodiumdiet=="low"])

summary(hypertension$bloodpressure[hypertension$sodiumdiet=="high"])
sd(hypertension$bloodpressure[hypertension$sodiumdiet=="high"])

library(tidyverse)
hypertension |> 
  group_by(sodiumdiet) |>
  summarize(
    Y_min = min(bloodpressure),
    Y_Q1 = quantile(bloodpressure, 0.25),
    Y_med = quantile(bloodpressure, 0.5),
    Y_Q3 = quantile(bloodpressure, 0.75),
    Y_max = max(bloodpressure),
    Y_IQR = Y_Q3 - Y_Q1,
    Y_mean = mean(bloodpressure),
    Y_sd = sd(bloodpressure)
  )

boxplot(hypertension$bloodpressure ~ hypertension$sodiumdiet, xlab="Sodium Diet", 
        ylab="Blood Pressure (mmHg)", main="Hypertension Experiment")

ggplot(hypertension, aes(x = sodiumdiet, y = bloodpressure)) +
  geom_boxplot() +
  ggtitle("Hypertension Experiment") +
  xlab("Sodium Diet") +
  ylab("Blood Pressure (mmHg)")

randomization.test <- function(response, treatment, Nsamps=10000, the.seed=500){
  ts <- rep(NA, Nsamps)
  n <- length(response)
  set.seed(the.seed)
  for(s in 1:Nsamps){
    permute <- sample.int(n)
    new.group <- treatment[permute]
    m1 <- mean(response[new.group==levels(treatment)[1]])
    m2 <- mean(response[new.group==levels(treatment)[2]])
    ts[s] <- (m1 - m2)
  }
  hist(ts, main="Randomization Distribution",
            xlab="Test Statistics (Difference in Mean Response)",
            ylab="Count")
  
  xbar1 <- mean(response[treatment==levels(treatment)[1]])
  xbar2 <- mean(response[treatment==levels(treatment)[2]])
  obsT <- (xbar1 - xbar2)
  p <- mean(abs(ts)>=abs(obsT))
  
  return(list(test_stats = ts, observed=obsT, p_value=p))
}
RT = randomization.test(hypertension$bloodpressure, hypertension$sodiumdiet)
RT
RT$test_stats[which(abs(RT$test_stats)>=abs(RT$observed))]
