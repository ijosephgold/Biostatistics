#2 independent groups, Wildtype and Mutant, are being compared for their chlorophyll content. Rows arent matched pairs, so an independent 2 sample comparison
#H_0: The means of Chlorophyll in the two groups are equal
#H_A: The means of Chlorophyll in the two groups are not equal
#First observe the data and calculate the means and sums, sample size, & standard error of the two groups. 
#Then, we will check the assumptions of normality and equal variance before deciding on the appropriate statistical test to use.

library(tidyverse) #for data wrangling
library(ggplot2) #for plotting
library(car) #for qqPlot function

Microalgaedata <- read.csv("data/Microalagae data.csv")
Microalgaedata

mean(Microalgaedata$Wildtype)
mean(Microalgaedata$Mutant, na.rm = TRUE) #na.rm mean remove missing values
sum(Microalgaedata$Wildtype)
sum(Microalgaedata$Mutant, na.rm = TRUE)

sample.size <- c(length(Microalgaedata$Wildtype), length(na.omit(Microalgaedata$Mutant)))
sample.size

standard.error <- c(sd(Microalgaedata$Wildtype)/sqrt(length(Microalgaedata$Wildtype)), sd(na.omit(Microalgaedata$Mutant))/sqrt(length(na.omit(Microalgaedata$Mutant))))
standard.error

qqPlot(Microalgaedata$Wildtype, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Wildtype Q-Q plot", col.lines = "red")
shapiro.test(Microalgaedata$Wildtype) #H_0: The data is normally distributed. H_A: The data is not normally distributed. If p-value is less than 0.05, the data is not normally distributed

qqPlot(na.omit(Microalgaedata$Mutant), distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Mutant Q-Q plot", col.lines = "red")
shapiro.test(na.omit(Microalgaedata$Mutant))

                      #data.frame creates table with 2 columns, strain {a label for each measurement} and chlorophyll {the actual measurement values}
Microalgaedata_long <- data.frame(
  strain = c(rep("Wildtype", length(Microalgaedata$Wildtype)), rep("Mutant", length(na.omit(Microalgaedata$Mutant)))),
  chlorophyll = c(Microalgaedata$Wildtype, na.omit(Microalgaedata$Mutant)))
Microalgaedata_long

#Bartlett test H0: The variances are equal
bartlett.test(chlorophyll ~ strain, data = Microalgaedata_long)

Microalgaedata_long$log_chlorophyll <- log(Microalgaedata_long$chlorophyll)
qqPlot(log_chlorophyll ~ strain, data = Microalgaedata_long, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Log-transformed Q-Q plot", col.lines = "red") # "==" means "is equal to?" true or false answers. == asks, <- assigns
shapiro.test(Microalgaedata_long$log_chlorophyl[Microalgaedata_long$strain == "Wildtype"]) #Square brackets tell R which entries to take out of the log_chlorophyll column to include only those in the shapiro test. In this case, we are taking out the log_chlorophyll values for the Wildtype strain.
shapiro.test(Microalgaedata_long$log_chlorophyll[Microalgaedata_long$strain == "Mutant"]) #
bartlett.test(log_chlorophyll ~ strain, data = Microalgaedata_long) #Bartlett test H0: The variances are equal
#No evidence the variances differ!

#H_0: The mean Log of Chlorophyll concentration in the two groups are equal
#H_A: The mean Log of Chlorophyll concentration in the two groups are not equal
t.test(log_chlorophyll ~ strain, data = Microalgaedata_long, var.equal = TRUE) 
#df = n_Wildtype + n_Mutant - 2 = 30 + 25 - 2 = 53
#p-value = 0.1193 > 0.05, so we fail to reject the null hypothesis.
#There is not enough evidence to conclude that the mean Log of Chlorophyll concentration is different between the two groups.