#Imports Spreadsheet
data <- read.csv("Data/Orange trees FA26 - Sheet1.csv")

#Opens Spreadsheet In New Tab
View(data)

#Calculates Mean Tree height for Orange Trees
meanheight <- data %>%
  summarise(mean(Tree.height..cm.))

#Calculates Standard Deviation for Tree height of Orange Trees
stdheight <- data %>%
  summarise(sd(Tree.height..cm.))

#Summary Function Above Generated Summary of Respective Data in a "mini-Table"

#Converts "mini-Table" Summary into Numeric Value
meanheight <- as.numeric(meanheight)

#calculate confidence limits
#first calculate summary stats for the entire data set
ybar <- mean(data$Tree.height..cm.)
SEMy <- sd(data$Tree.height..cm.)/sqrt(count(data))

#find the critical value of t with n-1 df
qt(p = 0.025, df = 131, lower.tail = FALSE)
lowerCL <- ybar-(SEMy*1.978239)
lowerCL
upperCL <- ybar+(SEMy*1.978239)
upperCL

#make a histogram w/ the confidence limits
hist(data$Tree.height..cm., breaks = 30, main = "", xlab= "Tree height cm", xlim = c(0,1600), ylim = c(0,25)) 
abline(v=lowerCL, col="blue")
abline(v=upperCL, col="blue")