# Load packages
library(car)
library(moments)

#Imports Spreadsheet
orangedata <- read.csv("Data/Orange trees FA26 - Sheet1.csv")

#Opens Spreadsheet In New Tab
View(orangedata)
orangedata #isolate collumns of interest: Circumference and Height of Trees
orangedata2 <- orangedata[,-c(1,2,4,5,6,7)]
orangedata2

#Calculate & renames Mean then Standard Deviation of Circumference for Orange Trees
meancir <- mean(orangedata2$Circumference..cm.)
meansd <- sd(orangedata2$Circumference..cm.)

n <- length(orangedata2$Circumference..cm.)
SEMmeancir <- meansd/sqrt(n) #Standard Error of the Mean: SE = s/sqrt(n), where s = sample standard deviation and n = sample size

#find the critical value of t with n-1 df
t_critical <- qt(1-0.05/2, df = n - 1) #: better & reusable way to find critical t value for 95% confidence interval
lowerCL <- meancir-(t_critical*SEMmeancir) #This calculates:[lower limit = ybar - t*SE]
lowerCL
upperCL <- meancir+(t_critical*SEMmeancir) #This calculates:[upper limit = ybar + t*SE]
upperCL

#Compare quantiles of the data to the quantiles of a normal distribution
qqPlot(orangedata2$Circumference, distribution = "norm", envelope = 0.95, main = "Normal Q-Q Plot of Tree Circumference")
skewness(orangedata2$Circumference..cm.) #near 0 means more symmetry, positive means longer right tail
kurtosis(orangedata2$Circumference..cm.) #normal distribution has kurtosis of 3 
shapiro.test(orangedata2$Circumference..cm.) #Shapiro p value greater than 0.5: no strong evidence against normality 

log_circumnference <- log(orangedata2$Circumference..cm.) #log transformation of circumference data
qqPlot(log_circumnference, distribution = "norm", envelope = 0.95, main = "Q-Q Plot of Log Circumference")
skewness(log_circumnference)
kurtosis(log_circumnference)
shapiro.test(log_circumnference)

sqrt_circumnference <- sqrt(orangedata2$Circumference..cm.) #square root transformation of circumference data
qqPlot(sqrt_circumnference, distribution = "norm", envelope = 0.95, main = "Q-Q Plot of Square Root Circumference")
skewness(sqrt_circumnference)
kurtosis(sqrt_circumnference)
shapiro.test(sqrt_circumnference)
#---------------------------------------------------------------
#Now for Height of Trees
height <- orangedata2$Tree.height..cm.
height <- height[!is.na(height)]

meanheight <- mean(height)
sdheight <- sd(height)
nheight <- sum(!is.na(height))

SEMheight <- sdheight/sqrt(nheight)
t_critical_height <- qt(0.975, df = nheight - 1)

lower_height <- meanheight - (t_critical_height * SEMheight)
lower_height
upper_height <- meanheight + (t_critical_height * SEMheight)
upper_height


qqPlot(height, distribution = "norm", envelope = 0.95, main = "Q-Q Plot of Tree Height", xlab = "Theoretical Normal Quantiles", ylab = "Observed Height")
skewness(height)
kurtosis(height)
shapiro.test(height)

log_height <- log(height) #log transformation of height data
qqPlot(log_height, distribution = "norm", envelope = 0.95, main = "Q-Q Plot of Log Height")
skewness(log_height)
kurtosis(log_height)
shapiro.test(log_height)

sqrt_height <- sqrt(height) #square root transformation of height data
qqPlot(sqrt_height, distribution = "norm", envelope = 0.95, main = "Q-Q Plot of Square Root Height")
skewness(sqrt_height)
kurtosis(sqrt_height)
shapiro.test(sqrt_height)

