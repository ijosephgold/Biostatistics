#Linear regression

#clear the environment
rm(list=ls())

#load libraries
library(tidyverse)
library(car) #for QQ-plot
library(ggplot2)

#let's analyze the species richness~biomass relationship

#load data
riffledata <- read.csv("Data/Stream data.csv")
riffledata

#first fit a linear model
model1 <- lm(NumberSpp~Biomass, data=riffledata, na.action="na.omit")
#plot model
plot(NumberSpp~Biomass, data=riffledata, col="blue", ylab="NumberSpp (mg)", xlab="Biomass (cm)")
abline(model1, col="blue") #adds fit from the model

#check the assumptions
#recall, we want to look at the residuals
#first test normality
model1res <- resid(model1)
qqp(model1res, "norm")
hist(model1res)
shapiro.test(model1res) #P < 0.05. It's a bit non-normal
#test for homogeneity of variance
#call the fitted y values
fitted(model1) #this is the predicted value of y, for each value of x
#plot the fitted (predicted) values against the residuals
plot(model1res~fitted(model1), ylab = "Residuals", xlab = "Predicted y") #it looks a bit cone shaped

#try a log transformation
riffledata$logBiomass<-log(riffledata$Biomass)
riffledata$logNumberSpp<-log(riffledata$NumberSpp)
riffledata

#create a new linear model with the log data
logmodel<-lm(logNumberSpp~logBiomass, data=riffledata, na.action="na.omit")

#check the assumptions of the log model
#test normality
logmodelres <- resid(logmodel)
qqp(logmodelres, "norm")
hist(logmodelres)
shapiro.test(logmodelres) #P > 0.05. This is better
#test for homogeneity of variance
# call the fitted y values
fitted(logmodel) 
#plot the fitted values against the residuals
plot(logmodelres~fitted(logmodel), ylab = "Residuals", xlab = "Predicted log y") #it looks better

#let's perfrom regression with the transformed data
summary(logmodel) #notice it also does a t-test (but t is calculated differently than in correlation)
#make an ANOVA table
anova(logmodel)

#make a scatterplot with trendline
ggplot(riffledata, aes(x=logBiomass, y=logNumberSpp)) +
  geom_point(size = 2) +
  xlab("Log (Biomass)") +
  ylab("Log (Species Richness)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
        geom_smooth(method = "lm") #adds regression line with 95% confidence intervals
