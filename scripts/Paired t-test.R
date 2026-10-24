#Paired t-test

#load libraries
library(tidyverse)
library(moments)
library(ggpubr)

#perform paired t-test on lab mice data

#create numeric data
#before drug
x <- as.numeric(c(98.5, 98.7, 98.8, 98.7, 98.6, 98.7, 98.7, 98.4, 98.5, 98.7,
                  98.2, 98.3, 98.5, 98.9, 98.2, 98.9, 98.8, 98.6, 98.4, 98.9, 
                  98.5, 98.5, 98.9, 98.4, 98.4, 98.7, 98.4, 98.7)) 
#after drug
y <- as.numeric(c(97.9, 98.0, 98.2, 98.3, 98.3, 98.4, 98.4, 98.5, 98.5, 98.5,
                  98.6, 98.6, 98.6, 98.7, 98.7, 98.7, 98.8, 98.9, 98.9, 99.0,
                  99.1, 99.2, 99.2, 99.3, 99.4, 99.5, 99.7, 99.9))

#combine samples into a data frame
drugdata <- data.frame(Before = x, After = y)
drugdata

#calculate differences
drugdata$Differences <- drugdata$Before - drugdata$After
drugdata #the above added a new column to the data frame

#calculate summary stats
mean <- mean(drugdata$Differences)
mean
stdev <- sd(drugdata$Differences)
stdev
n <- count(drugdata)
n
SEM <- stdev/sqrt(n) #SEM with n=28
SEM

#check model assumption of normality
hist(drugdata$Differences, breaks = 6, main = '') #looks good
hist(drugdata$Differences, breaks = 8, main = '') #bit of a negative skew
qqPlot(drugdata$Differences) #not too bad
skewness(drugdata$Differences) # >-1 and <1
kurtosis(drugdata$Differences) # <3
shapiro.test(drugdata$Differences) #P > 0.05

#calculate t
t <- (mean-0)/SEM
t

#perform test
#one tailed test (only concerned with temperatures elevating)
t.test(drugdata$Differences, alternative = c("less"), mu = 0, paired = FALSE)
#two-tailed test (temperatures could have gone in either direction)
t.test(drugdata$Differences, mu = 0, paired = FALSE)


#paired t-test on for chaparral abundance data

#load data
chaparraldata <- read.csv("Data/Chaparral data.csv")
chaparraldata
#to load data into a new window 
view(chaparraldata)

#check normality
hist(chaparraldata$Differences, main = '') #looks normal-ish
qqPlot(chaparraldata$Differences) #looks ok
skewness(chaparraldata$Differences) # <1
kurtosis(chaparraldata$Differences) # <3
shapiro.test(chaparraldata$Differences) #P > 0.05

#calculate summary stats
mean <- mean(chaparraldata$Differences)
mean
stdev <- sd(chaparraldata$Differences)
stdev
n <- count(chaparraldata)
n
SEM <- stdev/sqrt(n) #n=9
SEM

#calculate t
t <- (mean-0)/SEM #hypothesized mean difference of pairs = 0
t

#perform t-test
#two tailed test
t.test(chaparraldata$Desert, chaparraldata$Coast, paired = TRUE)
#one-tailed test (only tests if the abundances are higher near the coast)
t.test(chaparraldata$Desert,  alternative = c("less"), chaparraldata$Coast, paired = TRUE)

#create a paired box plot

#transform into long data
#the %>% operator is called a pipe - it tells R to take the data and then gather it
chaplong <- chaparraldata %>%
  gather(key = "site", value = "abundance", Desert, Coast)
chaplong #notice how the data is arranged now

#use ggpubr for plotting
ggpaired(chaplong, x = "site", y = "abundance", 
                order = c("Desert", "Coast"),
                ylab = "Relative abundance", xlab = "")

#or use boxplot() without lines
boxplot(chaparraldata$Coast, chaparraldata$Desert, names = c("Coast", "Desert"),
                ylab="Relative abundance")
