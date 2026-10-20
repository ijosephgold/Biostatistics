#H_0: Polarized Wings of female Heliconius butterflies had no efect on mate attraction 
#H_A: Polarized Wings of female Heliconius butterflies had an effect on mate attraction 
#Test of choice: T-test 
library(moments)
library(car)
library(ggplot2)
library(tidyverse)
library(ggpubr)

read.csv("Data/Butterfly wings.csv")
data <- read.csv("Data/Butterfly wings.csv")

#subtract Difference Depolarized - Polarized 
#Add collumn for Differences
data$Differences <- (data$Depolarized - data$Polarized)
data

hist(data$Differences, main = "Histogram of Differences in Mate Attraction", xlab = "Differences (Depolarized - Polarized)", col = "lightblue", border = "black")
qqPlot(data$Differences, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Q-Q plot of Differences in Mate Attraction", col.lines = "red")
skewness(data$Differences)  #if value for Skewness is 
kurtosis(data$Differences)
shapiro.test(data$Differences)

#transform into long data
#take the data then gather it
datalong <- data %>%
  gather(key = "Light", value = "Number of Visits", Polarized, Depolarized)
datalong
#useggPubr for plotting
ggpaired(datalong, x = "Light", y = "Number of Visits", 
         order = c("Polarized", "Depolarized"),
         ylab  = "Number of Visits", xlab = "Light treatment")

boxplot(data$Polarized, data$Depolarized, names = c("Polarized", "Depolarized"),
        col = "skyblue", xlab = "Light treatment", ylab = "Number of Visits")

t.test(data$Polarized, data$Depolarized, paired = TRUE)
#Based on the t-test results, we can reject the null hypothesis and conclude that there is a significant difference in mate attraction between polarized and depolarized wings.ox
