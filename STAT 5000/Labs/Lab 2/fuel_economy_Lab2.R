library(readr)
fuel <- read_csv("fuel_economy.csv", col_types = cols(Cylinders = col_factor(levels = c("4", "6"))))

library(tidyverse)
sum_stats = fuel |> 
  group_by(Cylinders) |>
  summarize(
    Y_n = n(),
    Y_mean = mean(Consumption.mpg),
    Y_sd = sd(Consumption.mpg)
  )
sum_stats

HT = t.test(Consumption.mpg~Cylinders, data=fuel, var.equal=T)
HT
