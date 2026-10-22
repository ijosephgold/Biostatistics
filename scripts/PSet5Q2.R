#ANOVA
#H_0: these conditions had no effect on plant size
#H_1: at least one condition had an effect on plant size

#load libraries
library(tidyverse)
library(car) #for QQ-plot
library(moments) #for skewness and kurtosis)
library(agricolae) #for Tukey's HSD
library(ggplot2)
library(lme4) #for random factors
library(vegan) #for perMANOVA

#Model I ANOVA
#test 

#load data
erodiumdata <- read.csv("Data/Erodium data.csv")
erodiumdata
inspect_df(erodiumdata) #check for missing data

#use data to calculate plant size as the product of height and width 
erodiumdata$productHW <- erodiumdata$Width*erodiumdata$Height
erodiumdata
inspect_df(erodiumdata)

#fit data to a linear model with Treatment = x and productHW = y
erodiumdata$Treatment <- as.factor(erodiumdata$Treatment)
lm <- lm(productHW~Treatment, data=erodiumdata) #takes the data as lm(y~x)
lm
anova(lm) #provides the ANOVA table
summary(lm) #provides extra details about the model

#check the assumptions
#recall, we want to look at the residuals to assess normality
lmres <- resid(lm)
qqPlot(lmres, "norm")
hist(lmres, xlab = "Log Product Height residuals", main = "Histogram for Log Transformed ")
boxplot(lmres, ylab = "Log Product Height residuals")
shapiro.test(lmres) #P > 0.05
skewness(lmres) 
kurtosis(lmres)
#next test equality of variances
bartlett.test(productHW ~ Treatment, data = erodiumdata) #P > 0.05

erodiumdata$logProduct <- log(erodiumdata$productHW) #log transform the data to see if it improves normality)
erodiumdata

#fit data to a linear model with Treatment = x and productHW = y
erodiumdata$Treatment <- as.factor(erodiumdata$Treatment)
lm <- lm(logProduct~Treatment, data=erodiumdata) #takes the data as lm(y~x)
lm

#check the assumptions
#recall, we want to look at the residuals to assess normality
lmres <- resid(lm)
qqPlot(lmres, "norm")
hist(lmres, xlab = "Log Product Height residuals", main = "")
boxplot(lmres, ylab = "Log Product Height residuals")
shapiro.test(lmres) #P > 0.05
skewness(lmres) 
kurtosis(lmres)
#next test equality of variances
bartlett.test(logProduct ~ Treatment, data = erodiumdata) #P > 0.05

anova(lm) #provides the ANOVA table
summary(lm) #provides extra details about the model
