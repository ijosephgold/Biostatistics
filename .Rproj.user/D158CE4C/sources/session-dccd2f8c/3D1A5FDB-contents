#Data transformations

#load libraries
library(car) #for QQ-plots
library(moments) #for analyzing skewness and kurtosis

#Let's look at the lab mouse temperature data to assess normality

#create numeric data
y <- as.numeric(c(97.9, 98.0, 98.2, 98.3, 98.3, 98.4, 98.4, 98.5, 98.5, 98.5,
                  98.6, 98.6, 98.6, 98.7, 98.7, 98.7, 98.8, 98.9, 98.9, 99.0,
                  99.1, 99.2, 99.2, 99.3, 99.4, 99.5, 99.7, 99.9))

#A histogram is one way to assess normality
hist(y, xlab = "Temperature", main = "") #looks pretty good
#What if we include more breaks?
hist(y, breaks = 8, xlab = "Temperature", main = "") #there's a bit of right skew
#We can also make a QQ plot
#The y-axis is the data's distribution and the x-axis is the expected normal distribution
#Look for a linear fit
qqnorm(y)
qqline(y) #most points are near the line
#Make a plot with the 95% confidence intervals
qqPlot(y, ylab = "Temperature") #most data are within the CLs 

#Analyze skewness and kurtosis 
skewness(y) #<1
kurtosis(y) #<3

#A Shapiro test is useful for assessing normality
#It tests the null hypothesis that the data is normally distributed
shapiro.test(y) #P = 0.698
#Seems like the data is normally distributed

#Look at the byssal threads data 
threads <- read.csv("Data/Byssal threads.csv")
hist(threads$Strength, xlab = "Strength", breaks = 8, main = "") #kind of right skewed
skewness(threads$Strength) #<1
kurtosis(threads$Strength) #<3
#Make a QQ-plot
qqPlot(threads$Strength, ylab = "Strength") #not great
#Try a Shapiro test
shapiro.test(threads$Strength) #P = 0.088
#It's safe to assume that normality is met, but can it be improved?

#Let's try a log transformation
threads$logstrength <- log(threads$Strength)
threads #notice logstrength is added to the data as a new column
hist(threads$logstrength, breaks = 8, xlab = "Log(Strength)",main = "") #kind of left skewed
shapiro.test(threads$logstrength) #P = 0.001
#That made it worse

#How about a square root transformation
threads$sqrtstrength <- sqrt(threads$Strength)
threads #sqrtstrength is another new column
hist(threads$sqrtstrength, xlab = "Sqrt(Strength)", breaks = 8, main = "") #not terrible
skewness(threads$sqrtstrength) #>-1
kurtosis(threads$sqrtstrength) #<3
qqPlot(threads$sqrtstrength, ylab = "Sqrt(Strength)")
shapiro.test(threads$sqrtstrength) #P = 0.22
#This is a slight improvement

#Look at kelp Bass gonad mass data
gonads <- read.csv("Data/Gonad data.csv")
hist(gonads$gonad_mass, xlab = "Gonad mass", breaks = 20, main = "") #a lot of right skew
skewness(gonads$gonad_mass) #>1
kurtosis(gonads$gonad_mass) #>3
qqPlot(gonads$gonad_mass, ylab = "Gonad mass")
shapiro.test((gonads$gonad_mass)) #P<<0.0001
#The data is far from normally distributed

#Try a log transformation
gonads$logmass <- log(gonads$gonad_mass)
hist(gonads$logmass, xlab = "Log(Gonad mass)", breaks = 20) #better but its left skewed
skewness(gonads$logmass) #>-1
kurtosis(gonads$logmass) #>3
qqPlot(gonads$logmass, ylab = "Log(Gonad mass)")
shapiro.test(gonads$logmass) #P<<0.0001
#A slight improvement

#Square root transformation
gonads$sqrtmass <- sqrt(gonads$gonad_mass)
hist(gonads$sqrtmass, xlab = "Sqrt(Gonad mass)", breaks = 20) #kind of better 
skewness(gonads$sqrtmass) #<1
kurtosis(gonads$sqrtmass) #<3
qqPlot(gonads$sqrtmass, ylab = "Sqrt(Gonad mass)")
shapiro.test(gonads$sqrtmass) #P<<0.0001
#Better than the raw data and the log transformation
