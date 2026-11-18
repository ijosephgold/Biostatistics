```r
load_biostats_packages()
#load libraries
library(tidyverse)
library(car)
library(ggplot2)
library(GGally)
library(psych)
library(moments)
#H_0: Beak width and beak depth have no linear association in the population; the correlation is zero (Ï = 0)
#H_A: Beak width and beak depth have a linear association; the correlation is not zero (Ï â‰  0)
penguins <- read.csv("data/penguin data - Sheet1.csv")
View(penguins)

#penguins$Beak.width <- as.numeric(penguins$Beak.width)
qqPlot(penguins$Beak.width, "norm",ylab = "", main = "Q-Q plot of Beak width")
skewness(penguins$Beak.width) #near 0 = more symmetrical, + = longer right-tail, - = longer left-tail 
kurtosis(penguins$Beak.width) #normal distribution ~ 3, >3 = more peaked, <3 = flatter
shapiro.test(penguins$Beak.width) #p-value > 0.05 = no strong evidence against normality, p-value <= 0.05 = evidence against normality
hist(penguins$Beak.width, main = "Histogram of Beak width", xlab = "Beak width (mm)", ylab = "Frequency", col = "lightblue", border = "black")

#penguins$Beak.width <- as.numeric(penguins$Beak.width)
qqPlot(penguins$Beak.depth, "norm", ylab = "", main = "Q-Q plot of Beak depth")
skewness(penguins$Beak.depth) #near 0 = more symmetrical, + = longer right-tail, - = longer left-tail 
kurtosis(penguins$Beak.depth) #normal distribution ~ 3, >3 = more peaked, <3 = flatter
shapiro.test(penguins$Beak.depth) #p-value > 0.05 = no strong evidence against normality, p-value <= 0.05 = evidence against normality
hist(penguins$Beak.depth, main = "Histogram of Beak depth", xlab = "Beak depth (mm)", ylab = "Frequency", col = "lightblue", border = "black")

#scatterplot the relationship between beak width and beak depth
plot(penguins$Beak.width ~ penguins$Beak.depth, data = penguins,
     xlab = "Beak depth (mm)",
     ylab = "Beak width (mm)",
     main = "Beak width by Beak depth in African penguins",
     pch = 19)
abline(lm(penguins$Beak.width ~ penguins$Beak.depth, data = penguins), col = "red")

cor.test(penguins$Beak.depth, penguins$Beak.width, method = "pearson")

#penguins_long <- penguins %>%
#  gather(key = "Trait", value = "Size", Beak.width, Beak.depth)
#penguins_long 
#view(penguins_long)
#bartlett.test(Size ~ Trait, data = penguins_long)


#plot(penguins$Beak.width, penguins$Beak.depth,
#     xlab = "Beak width (mm)",
#     ylab = "Beak depth (mm)",
#     main = "Beak width and depth in African penguins",
#     pch = 19)
#abline(lm(penguins$Beak.depth ~ penguins$Beak.width), col = "blue", lwd = 2)

#model <- lm(penguins$Beak.depth ~ penguins$Beak.width)

#plot(model, which = 1)  
#Plot 1: Residuals vs. Fitted. residual is (depth y_i) - (depth y-hat), looking for residuals scattered random around 0, curve ~suggest relationship isnâ€™t linear; funnel shape ~suggest residual spread changes across fitted values
#plot(model, which = 2)  
#Plot 2: Normal Q-Q. compares residuals with what weâ€™d expect if they were normally distributed. Points roughly following line suggest normality is reasonable, visually complements Shapiroâ€“Wilk test ran

#shapiro.test(resid(model))


#penguins$differences <- penguins$Beak.depth - penguins$Beak.width
#var(penguins$differences)
#qqPlot(penguins$differences)
#shapiro.test(penguins$differences)

```
