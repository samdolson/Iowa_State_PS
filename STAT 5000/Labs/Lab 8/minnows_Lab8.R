library(readr)
minnows <- read_table("minnows.txt", 
                      col_names = c("copper", "zinc", "protein"), 
                      col_types = cols(copper = col_factor(levels = c("0","150")), 
                              zinc = col_factor(levels = c("0", "750", "1500"))))

options(contrasts = c("contr.sum", "contr.sum"))
sumto0.anova <- aov(protein ~ copper + zinc + copper*zinc, data = minnows)
summary(sumto0.anova)
sumto0.anova$coefficients
model.matrix(sumto0.anova)

options(contrasts = c("contr.treatment", "contr.treatment"))
baseline.anova <- aov(protein ~ copper + zinc + copper*zinc, data = minnows)
summary(baseline.anova)
baseline.anova$coefficients
model.matrix(baseline.anova)

options(contrasts = c("contr.sum", "contr.sum"))
summary(lm(protein ~ copper + zinc + copper*zinc, data = minnows))

interaction.plot(minnows$copper, minnows$zinc, minnows$protein,
                 main="Minnows Experiment", 
                 xlab="Zinc Concentration",
                 ylab="Mean Protein Level",
                 trace.label="Copper\n Conc.")
interaction.plot(minnows$zinc, minnows$copper, minnows$protein,
                 main="Minnows Experiment", 
                 xlab="Copper Concentration",
                 ylab="Mean Protein Level",
                 trace.label="Zinc\n Conc.")

library(emmeans)
simple.effects <- emmeans(sumto0.anova, c("copper", "zinc"))
pairs(simple.effects, adjust=NULL)

copper.effects <- emmeans(sumto0.anova, "copper")
pairs(copper.effects, adjust="tukey")

zinc.effects <- emmeans(sumto0.anova, "zinc")
pairs(zinc.effects, adjust="tukey")

plot(sumto0.anova$fitted.values, sumto0.anova$residuals,
     main="Diagnostic Plot for Minnows Experiment",
     xlab="Fitted Values", ylab="Residuals")
plot(as.numeric(minnows$copper), sumto0.anova$residuals,
     main="Diagnostic Plot for Minnows Experiment",
     xlab="Copper Level", ylab="Residuals")
plot(as.numeric(minnows$zinc), sumto0.anova$residuals,
     main="Diagnostic Plot for Minnows Experiment",
     xlab="Zinc Level", ylab="Residuals")

qqnorm(sumto0.anova$residuals)
qqline(sumto0.anova$residuals, col="red")
shapiro.test(sumto0.anova$residuals)
library(moments)
mean(sumto0.anova$residuals)
median(sumto0.anova$residuals)
skewness(sumto0.anova$residuals)
kurtosis(sumto0.anova$residuals)-3



