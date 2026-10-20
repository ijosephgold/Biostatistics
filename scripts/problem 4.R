#2 independent groups, Wildtype and Mutant, are being compared for their chlorophyll content. Rows arent matched pairs, so an independent 2 sample comparison
#H_0: The means of Chlorophyll in the two groups are equal
#H_A: The means of Chlorophyll in the two groups are not equal
#First observe the data and calculate the means and sums, sample size, & standard error of the two groups. 
#Then, we will check the assumptions of normality and equal variance before deciding on the appropriate statistical test to use.

library(tidyverse) #for data wrangling
library(ggplot2) #for plotting
library(car) #for qqPlot function
library(moments)

Microalgaedata <- read.csv("data/Microalagae data.csv")
Microalgaedata
inspect_df(Microalgaedata)

Microalgaedata_long <- Microalgaedata %>%
  gather(key = "Genotype", value = "Concentration", Wildtype, Mutant, na.rm = TRUE)
Microalgaedata_long
inspect_df(Microalgaedata_long)

qqPlot(Microalgaedata_long$Concentration, ylab = "Concentration", distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Concentration Q-Q plot", col.lines = "red")
shapiro.test(Microalgaedata_long$Concentration) #H_0: The data is normally distributed. H_A: The data is not normally distributed. If p-value is less than 0.05, the data is not normally distributed
bartlett.test(Concentration ~ Genotype, data = Microalgaedata_long)

Microalgaedata_long$log_Concentration <- log(Microalgaedata_long$Concentration)
Microalgaedata_long
qqPlot(Microalgaedata_long$log_Concentration, ylab = "Log Concentration", distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Log Concentration Q-Q plot", col.lines = "red")
shapiro.test(Microalgaedata_long$log_Concentration)  #H_0: The data is normally distributed. H_A: The data is not normally distributed. If p-value is less than 0.05, the data is not normally distributed
bartlett.test(log_Concentration ~ Genotype, data = Microalgaedata_long)

hist(Microalgaedata_long$Concentration, main = "Histogram of Concentration", xlab = "Concentration", col = "lightblue", border = "black")
hist(Microalgaedata_long$log_Concentration, main = "Histogram of Log Concentration", xlab = "Log Concentration", col = "lightblue", border = "black")
boxplot(Microalgaedata$Wildtype, Microalgaedata$Mutant, names = c("Wildtype", "Mutant"),
        col = "skyblue", xlab = "Strain", ylab = "Concentration")
boxplot(log_Concentration ~ Genotype, data = Microalgaedata_long, col = "skyblue", xlab = "Strain", ylab = "Log Concentration")
t.test(log_Concentration ~ Genotype, data = Microalgaedata_long, var.equal = TRUE) 
