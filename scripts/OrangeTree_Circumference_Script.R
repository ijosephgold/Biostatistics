#Imports Spreadsheet
data <- read.csv("Data/Orange trees FA26 - Sheet1.csv")

#Opens Spreadsheet In New Tab
View(data)

#Calculates Mean Circumference for Orange Trees
meancir <- data %>%
  summarise(mean(Circumference..cm.))

#Calculates Standard Deviation for Circumference of Orange Trees
stdcir <- data %>%
  summarise(sd(Circumference..cm.))

#Summary Function Above Generated Summary of Respective Data in a "mini-Table"

#Converts "mini-Table" Summary into Numeric Value
meancir <- as.numeric(meancir)

#calculate confidence limits
#first calculate summary stats for the entire data set
ybar <- mean(data$Circumference..cm.)
SEMy <- sd(data$Circumference..cm.)/sqrt(count(data))

#find the critical value of t with n-1 df
qt(p = 0.025, df = 131, lower.tail = FALSE)
lowerCL <- ybar-(SEMy*1.978239)
lowerCL
upperCL <- ybar+(SEMy*1.978239)
upperCL

#make a histogram w/ the confidence limits
hist(data$Circumference..cm., breaks = 30, main = "", xlab= "Circumference", ylim = c(0,25)) 
abline(v=lowerCL, col="blue")
abline(v=upperCL, col="blue")