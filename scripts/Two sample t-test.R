#Two-sample t-test

#load libraries
library(tidyverse)
library(car)
library(moments)
library(ggplot2)

#test barnacle recruitment in open ocean vs. bays

#load data
barnacledata <- read.csv("Data/Barnacle data.csv")
barnacledata

#look at variance
#openocean
openvar <- barnacledata %>%
  filter(Site == "openocean") %>%
  summarize(var=var(recruitment, na.rm=TRUE))
openvar
#bays
bayvar <- barnacledata %>%
  filter(Site == "bay") %>%
  summarize(var=var(recruitment, na.rm=TRUE))
bayvar
#they're not the same, but within the same order of magnitude

#try a Bartlett test
bartlett.test(recruitment~Site, data=barnacledata) #P > 0.05

#check normality
#make a histogram
hist(barnacledata$recruitment, main= '') #very right skewed
shapiro.test(barnacledata$recruitment) #P << 0.05
skewness(barnacledata$recruitment) #1.01
kurtosis(barnacledata$recruitment) #4.10
qqPlot(barnacledata$recruitment)
#data is not normally distributed

#try a sqrt transformation (the data has zeros)
barnacledata$sqrtdata <- sqrt(barnacledata$recruitment)
barnacledata

#check normality
#make histogram
hist(barnacledata$sqrtdata, main = '') #very right skewed
shapiro.test(barnacledata$sqrtdata) #P = 0.003, this is better but still < 0.05
skewness(barnacledata$sqrtdata) #good
kurtosis(barnacledata$sqrtdata) #good
qqPlot(barnacledata$sqrtdata)
#data is still not normally distributed but the transformation improved it

#try a log transformation
#all the data get a +1 since there are zeros
barnacledata$plusone <- (barnacledata$recruitment)+1
barnacledata
barnacledata$logdata <- log(barnacledata$plusone)
barnacledata

#check normality
#make histogram
hist(barnacledata$logdata) #very weird
shapiro.test(barnacledata$logdata) #P << 0.05, this is worse
#not even gonna bother looking at skewness & kurtosis
qqPlot(barnacledata$logdata)
#log transformation doesn't help at all

#run two sample t-tests
#for a two sample test, the groups have to be specified as independetvar~categoricalvar
#with sqrt transformed data
t.test(sqrtdata~Site, var.equal=TRUE, data=barnacledata)
#with regular data
t.test(recruitment~Site, var.equal=TRUE, data=barnacledata)
#notice that these tests assume equal variance "var.equal = TRUE"
#if they weren't, instead run Welch's test by setting "var.equal = FALSE"
t.test(sqrtdata~Site, var.equal=FALSE, data=barnacledata)
#notice that the df and p-value are not the same

#try a Mann-Whitney U-test with raw data
wilcox.test(recruitment~Site, data=barnacledata)

#make a barplot with means +/- standard error bars
#first put together summary stats
graphdata <- barnacledata %>%
  group_by(Site) %>%
  summarize(meanrecruitment=mean(recruitment, na.rm=TRUE), se=sd(recruitment, na.rm=TRUE)/sqrt(30))
graphdata

#ggplot is best for making a plot with error bars

#first look at how to do this step-by-step
#notice that the object p (the plot) can be built upon with p+'s

#define what x and y variables to use to make the graph
p<-ggplot(graphdata, aes(x=Site, y=meanrecruitment)) 
#gts rid of the gray background and the useless gridlines
p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))
#specify what colors should the bars be. Why not think Dodger blue?
p<-p + geom_bar(stat="identity", color="dodgerblue2", fill="dodgerblue2", position="dodge", size=0.6) 
#label the axes
p<-p + labs(x="Site", y="Recruitment", size = 1) 
#add error bars
p<-p + geom_errorbar(aes(ymax=meanrecruitment+se, ymin=meanrecruitment-se), position=position_dodge(0.9), width=0.1) 
#print the final product
print(p)

#it's possible to do all of the above by instead adding the parts together step-wise
ggplot(graphdata, aes(x=Site, y=meanrecruitment)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_bar(stat="identity", color="dodgerblue2", fill="dodgerblue2", position="dodge", size=0.6) +
  labs(x="Site", y="Recruitment", size = 1) +
  geom_errorbar(aes(ymax=meanrecruitment+se, ymin=meanrecruitment-se), position=position_dodge(0.9), width=0.1) 
