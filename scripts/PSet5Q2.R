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
hist(lmres, xlab = "Product Height residuals", main = "Histogram")
boxplot(lmres, ylab = "Product Height residuals")
shapiro.test(lmres) #P > 0.05
skewness(lmres) 
kurtosis(lmres)
#next test equality of variances
bartlett.test(productHW ~ Treatment, data = erodiumdata) #P > 0.05
#hist(resid(lm))


#log transform the data to see if it improves normality
erodiumdata$logProduct <- log(erodiumdata$productHW) #log transform the data to see if it improves normality)
erodiumdata

#fit data to a linear model with Treatment = x and productHW = y
erodiumdata$Treatment <- as.factor(erodiumdata$Treatment)

loglm <- lm(logProduct~Treatment, data=erodiumdata) #takes the data as lm(y~x)
loglm

#check the assumptions
#recall, we want to look at the residuals to assess normality
loglmres <- resid(loglm)

qqPlot(loglmres, "norm")
hist(loglmres, xlab = "Log Product Height residuals", main = "")
boxplot(loglmres, ylab = "Log Product Height residuals")
shapiro.test(loglmres) #P > 0.05
skewness(loglmres) 
kurtosis(loglmres)
#next test equality of variances
bartlett.test(logProduct ~ Treatment, data = erodiumdata) #P > 0.05

anova(loglm) #provides the ANOVA table
summary(loglm) #provides extra details about the model


#prepare data for plotting
graphdata <- erodiumdata %>%
  group_by(Treatment) %>%
  summarize(mean=mean(productHW, na.rm=TRUE), se=sd(productHW, na.rm=TRUE)/sqrt(length(na.omit(productHW))))
graphdata
#rearrange order of the levels on the x-axis
graphdata$Treatment <- factor(graphdata$Treatment, levels = c("control", "shade", "wind", "3XH2o"))

#make a barplot with ggplot()
ggplot(data=graphdata, aes(x=Treatment, y=mean, fill=Treatment)) + 
  theme_bw() + 
  theme(axis.text.x=element_text(color="black", size=14), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=16, face="bold"), 
        axis.title.y = element_text(color="black", size=16, face="bold"),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_bar(colour="black", fill="grey", width=0.8, stat="identity") + 
  guides(fill=FALSE) + 
  ylab("Product Height x Width") +
  xlab("Treatment Type") +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), stat="identity", 
                position=position_dodge(width=0.9), width=0.1) +
  geom_text(aes(label=c("a", "bc", "b", "c"), y = mean + 1))
#note that the labels are assigned in the order of levels in "graphdata"
