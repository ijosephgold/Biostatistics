#Contingency tables
#Imports Spreadsheet
data <- read.csv("Data/Ladybird data.csv")

#Opens Spreadsheet In New Tab
View(data)

#load libraries 
library(tidyverse)
library(DescTools) #for G test
library(ggplot2) #for plotting

#analyze dispersal and reproduction data

#create data
redrural <- 70
redindustrial <- 85
blackrural <- 30
blackindustrial <- 115

# Make contingency table using ONLY the observed values
contable_obs <- matrix(c(redrural, blackrural,
                         redindustrial, blackindustrial),
                       nrow = 2, byrow = TRUE)

# Name rows and columns
rownames(contable_obs) <- c("Rural", "Industrial")
colnames(contable_obs) <- c("Red", "Black")

addmargins(contable_obs)

# Calculate expected values
contable_exp <- outer(rowSums(contable_obs),
                      colSums(contable_obs)) / sum(contable_obs)

contable_exp

# Perform G-test
GTest(contable_obs)

#make plot
names <- matrix(c("Rural", "Ind"))
barplot(contable_obs, beside = TRUE, legend.text = rownames(contable_obs), xlab = "Color", ylab = "Counts", args.legend = list(x = "topright"))