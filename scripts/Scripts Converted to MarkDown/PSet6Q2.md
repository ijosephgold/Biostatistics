```r
#clear the environment
rm(list=ls())

#load libraries
library(tidyverse)
library(car)
library(ggplot2)
library(GGally)
library(psych)

edata <- read.csv("data/Environmental data.csv")
view(edata)
str(edata)

data_subset <-  edata[, c("PctCover", "Elevation", "TPI", "TWI", "MaxTemp")]

#create the scatter plot matrix
pairs(data_subset, 
      main = "SPLOM",
      pch = 20,                                     
      cex = 1.5,                                    
      col = "#69b3a2"                               
)
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
cor.test(edata$PctCover, edata$Elevation, method="pearson", na.rm=TRUE)
cor.test(edata$PctCover, edata$TPI, method="pearson", na.rm=TRUE)
cor.test(edata$PctCover, edata$TWI, method="pearson", na.rm=TRUE)
cor.test(edata$PctCover, edata$MaxTemp, method="pearson", na.rm=TRUE)
cor.test(edata$Elevation, edata$TPI, method="pearson", na.rm=TRUE)
cor.test(edata$Elevation, edata$TWI, method="pearson", na.rm=TRUE)
cor.test(edata$Elevation, edata$MaxTemp, method="pearson", na.rm=TRUE) #***
cor.test(edata$TPI, edata$TWI, method="pearson", na.rm=TRUE)
cor.test(edata$TPI, edata$MaxTemp, method="pearson", na.rm=TRUE)
cor.test(edata$TWI, edata$MaxTemp, method="pearson", na.rm=TRUE)


```
