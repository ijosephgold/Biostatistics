#Correlation

#clear the environment
rm(list=ls())

#load libraries
library(tidyverse)
library(car)
library(ggplot2)
library(GGally)
library(psych)

#analyze desert ant nest distributions 
#test for correlation between mean distance to vegetation and neighbor distances

#load data
antdata <- read.csv("Data/Ant data.csv")
View(antdata)
str(antdata) #look at structure of variables in data frame

#there's alot of info here
#let's work with the columns "Mean" and "Neighbors"

#check assumptions of normality
#start with mean data
hist(antdata$Mean)
shapiro.test(antdata$Mean)
qqPlot(antdata$Mean)
#the data for mean distance is not normally distributed

#try a log transformation
antdata$logmean <- log(antdata$Mean)
hist(antdata$logmean)
shapiro.test(antdata$logmean)
qqPlot(antdata$logmean)
#log transformation improves normality

#what about square root transformation
antdata$sqrtmean <- sqrt(antdata$Mean)
hist(antdata$sqrtmean)
shapiro.test(antdata$sqrtmean)
qqPlot(antdata$sqrtmean)
#not as good as the log transformation

#how about neighbor distance data
hist(antdata$Neighbors)
shapiro.test(antdata$Neighbors) 
qqPlot(antdata$Neighbors)
#also pretty bad

#let's try a log transformation
antdata$logneighbors <- log(antdata$Neighbors)
hist(antdata$logneighbors)
shapiro.test(antdata$logneighbors) 
qqPlot(antdata$logneighbors)
#log transformation improves normality

#what about square root transformation
antdata$sqrtneighbors <- sqrt(antdata$Neighbors)
hist(antdata$sqrtneighbors)
shapiro.test(antdata$sqrtneighbors)
qqPlot(antdata$sqrtneighbors)
#not as good as the log transformation

#let's try Pearson's correlation with and without the log transformed data
cor.test(antdata$logneighbors, antdata$logmean, method="pearson", na.rm=TRUE)
cor.test(antdata$Neighbors, antdata$Mean, method="pearson", na.rm=TRUE)
#it's a highly significant association either way
#the log transformation has slightly better model fit (r = 0.80 vs. 0.77)

#how about Spearman's rank test (non-parametric) with raw data
cor.test(antdata$Neighbors, antdata$Mean, method="spearman", na.rm=TRUE)
#it gives a warning about ties (the ranks don't fit perfectly)
#interpret the P value as an approximation

#make a scatterplot
ggplot(antdata, aes(x=logmean, y=logneighbors)) +
  geom_point(size = 2) +
  xlab("Log (Mean distance to vegetation)") +
  ylab("Log (Mean distance to neighbors)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
#we could color code the points by site
ggplot(antdata, aes(x=logmean, y=logneighbors, col = Site)) +
  geom_point(size = 2) +
  xlab("Log (Mean distance to vegetation)") +
  ylab("Log (Mean distance to neighbors)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())

#let's make a SPLOM with Transverse range transect data
#we'll look at shrub diversity, evenness, percent cover, and species richness

#load data
transects <- read.csv("Data/Transverse transect data.csv")
View(transects)
str(transects) #look at structure of variables in data frame

#we should probably check the model assumptions
#let's skip that here as the code to do so is available above

#select columns
data_subset <- transects[, c("Diversity", "Evenness", "Richness", "Cover")]

#create the scatter plot matrix
pairs(data_subset, 
      main = "SPLOM",
      pch = 20,                                     
      cex = 1.5,                                    
      col = "#69b3a2"                               
)

#alternate method using ggpairs() from library(GGgally) that looks a little better
#the diagonal plots are histograms
ggpairs(data_subset, 
        title = "")

#yet another method using pair.panels() from library(psych)
#this gives histograms and ellipses
pairs.panels(data_subset,
             method = "pearson", # Correlation method (e.g., pearson, spearman)
             hist.col = "#00AFBB", # Color for histograms
             density = TRUE,      # Add density plots
             ellipses = TRUE      # Add correlation ellipses
)

#create the individual Pearson's correlation tests
#notice that the results match those shown above
cor.test(transects$Diversity, transects$Evenness, method="pearson", na.rm=TRUE)
cor.test(transects$Diversity, transects$Richness, method="pearson", na.rm=TRUE)
cor.test(transects$Diversity, transects$Cover, method="pearson", na.rm=TRUE)
cor.test(transects$Evenness, transects$Richness, method="pearson", na.rm=TRUE)
cor.test(transects$Evenness, transects$Cover, method="pearson", na.rm=TRUE)
cor.test(transects$Cover, transects$Richness, method="pearson", na.rm=TRUE)
