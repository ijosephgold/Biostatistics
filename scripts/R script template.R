#R script template

#Getting a new project going
#It's imperative to have an organized directory to make R work smoothly

#FIRST
#Create three new folders within your project's directory
#This can be done outside of R studio, or in the Files window at lower right
#Name the folders:
#Data
#Script
#Output

#You can copy the first two scripts into the Script folder outside of R studio

#The first line of an R Script should be some kind of title
#Followed by the author's name and date
#For example:

#New R Script
#Created by you
#Date

#Next provide a brief description of what the code is going to do
#We'll work with the Mytilus byssal thread break force data here for practice

#Often times, a script will begin by sweeping the environment clean
rm(list=ls())
#This can also be done by clicking the broom above the environment window

#SECOND
#Load libraries
#Libraries are software packages
#Click on the Packages tab at lower right to see a list of them
#You can install them by clicking on the box to the left
#Or by using the function install.packages()
#For example:
install.packages('tidyverse')
#Notice that the package is being installed in the console
#You can usually ignore this until it's done 
#Now load the package with the library() function
library(tidyverse)
#The tidyverse package has many useful functions for data cleaning
#We'll use lots of other packages eventually
#For now the tidyverse will help us bring an Excel spreadsheet into R

#THIRD
#Load data
#Let's start with loading the spreadsheet created in Excel
#We are going to make it a new object
data <- read.csv("Data/Byssal threads(sheet1).csv")
#read.csv calls the spreadsheet
#The directory must be specified exactly or it won't work
#You can find your working directory with getwd()
getwd()
#View the data in a new tab with View()
View(data)
#Or just run it in the console
data
#Look at the first few rows only
head(data)
#Look at the data's structure
str(data)

#FOURTH
#We can calculate summary statistics by site
#We need the %>% from the tidyverse (this is called a pipe)
#Pristine site
prmeanstr <- data %>%
  filter(Site=="Pristine") %>%
  summarise(mean(Strength)) #returns the mean strength at the pristine site
prstdstr <- data %>%
  filter(Site=="Pristine") %>%
  summarise(sd(Strength)) #returns the stdev strength at the pristine site
#Polluted site
plmeanstr <- data %>%
  filter(Site=="Polluted") %>%
  summarise(mean(Strength)) #returns the mean strength at the polluted site
plstdstr <- data %>%
  filter(Site=="Polluted") %>%
  summarise(sd(Strength)) #returns the stdev strength at the polluted site

#We can make a simple graph
#first set means as numeric
prmeanstr <- as.numeric(prmeanstr)
plmeanstr <- as.numeric(plmeanstr)
#concatenate them into a variable y
y <- c(prmeanstr, plmeanstr)
#name the levels of x
names <- matrix(c("Pristine", "Polluted"))
#make a data frame
graphdata <- data.frame(Site = names, Strength = y)
#use barplot() for graphing
barplot(height = graphdata$Strength, names.arg = graphdata$Site, 
        ylim = c(0,1000), xlab = "Site", ylab = "Mean Strength (g)")

#Add error bars with the standard deviations using the arrows() function
#first make stdev numeric
prstdstr <- as.numeric(prstdstr)
plstdstr <- as.numeric(plstdstr)
#concatenate them into an object
stdev <- c(prstdstr, plstdstr)
#redefine the plot with greater y axis limit that fits bars
barploter <- barplot(height = graphdata$Strength, names.arg = graphdata$Site, 
                     ylim = c(0,1500), xlab = "Site", ylab = "Mean Strength (g)")
#add bars
arrows(x0 = barploter, y0 = y - stdev,
       x1 = barploter, y1 = y + stdev,
       angle = 90, code = 3, length = 0.05)

#calculate confidence limits
#first calculate summary stats for the entire data set
ybar <- mean(data$Strength)
SEMy <- sd(data$Strength)/sqrt(count(data))
#find the critical value of t with n-1 df
qt(p = 0.05, df = 24, lower.tail = FALSE)
lowerCL <- ybar-(SEMy*1.710882)
lowerCL
upperCL <- ybar+(SEMy*1.710882)
upperCL

#make a histogram w/ the confidence limits
hist(data$Strength, breaks = 8, main = "", xlab= "Strength")
abline(v=lowerCL, col="blue")
abline(v=upperCL, col="blue")
