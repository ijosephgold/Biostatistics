#One sample t-test

#load tidyverse
library(tidyverse)

#perform one sample t-test on lab mice data

#create numeric data
y <- as.numeric(c(97.9, 98.0, 98.2, 98.3, 98.3, 98.4, 98.4, 98.5, 98.5, 98.5,
                  98.6, 98.6, 98.6, 98.7, 98.7, 98.7, 98.8, 98.9, 98.9, 99.0,
                  99.1, 99.2, 99.2, 99.3, 99.4, 99.5, 99.7, 99.9))

#make a histogram to see the spread of the data
hist(y, main = "")
#change the number of bars with breaks
hist(y, breaks = 8, main = "")

#a boxplot shows the data in a different way
boxplot(y, ylim = c(97.8,100.1), ylab ="Temperature (F)") 
# I used ylim() to set the y axis limits below and above 0.1 of the most extreme values

#calculate summary stats
ybar <- mean(y) #average
ybar
stdevy <- sd(y) #standard deviation
stdevy
ydata <- as.data.frame(y) #set up a data frame to check the count
n <- count(ydata)
n
SEMy <- stdevy/sqrt(n) #SEM with n=28
SEMy

#calculate t
t <- (ybar-98.6)/SEMy
t

#perform test
t.test(y, alternative = c("greater"), mu = 98.6, paired = FALSE)

#calculate confidence limits
#one tailed
#critical value of t with df=27 is 1.703
lowerbound <- ybar-(SEMy*1.703)
lowerbound
#two tailed
#critical value of t with  df=27 is 2.052
lowerCL <- ybar-(SEMy*2.052)
lowerCL
upperCL <- ybar+(SEMy*2.052)
upperCL

#make a histogram w/ the confidence limits
#one tail
hist(y, breaks = 8, xlab = "Temperature", main = "")
abline(v=lowerbound, col="blue")
#two tails
hist(y, breaks = 8, xlab = "Temperature", main = "")
abline(v=lowerCL, col="blue")
abline(v=upperCL, col="blue")

#or a boxplot with notches showing two-tailed CL's
boxplot(y, ylim = c(97.8,100.1), notch = TRUE, ylab ="Temperature (F)") 
