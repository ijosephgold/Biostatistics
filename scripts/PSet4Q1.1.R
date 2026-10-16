#Imports Spreadsheet
orangedata <- read.csv("../Data/Orange trees FA26 - Sheet1.csv")

#Opens Spreadsheet In New Tab
View(orangedata)
orangedata
#orangedata2 <- orangedata[,-c(1,2,4,5,6,7)])

#orangedata2

#Calculates Mean Circumference for Orange Trees
meancir <- orangedata %>%
  summarise(mean(Circumference..cm.))
#Alt cleaner version: meancir <- mean(data$Circumference..cm., na.rm = TRUE) (na.rm = TRUE ignores missing values)
#Using Alt cleaner version negates the need to convert "mini-Table" Summary into Numeric Value below

#Calculates Standard Deviation for Circumference of Orange Trees
stdcir <- orangedata %>%
  summarise(sd(Circumference..cm.))
#Alt cleaner version: stdcir <- sd(orangedata$Circumference..cm., na.rm = TRUE)

#Summary Function Above Generated Summary of Respective Data in a "mini-Table"

#Converts "mini-Table" Summary into Numeric Value
meancir <- as.numeric(meancir)

#calculate confidence limits
#first calculate summary stats for the entire data set
ybar <- mean(orangedata$Circumference..cm.) #sample mean

#n <- sum(!is.na(orangedata$Circumference..cm.))
#SEMy <- sd(orangedata$Circumference..cm., na.rm = TRUE) / sqrt(n)
SEMy <- sd(orangedata$Circumference..cm.)/sqrt(count(orangedata)) #Standard Error of the Mean: SE = s/sqrt(n), where s = sample standard deviation and n = sample size

#find the critical value of t with n-1 df
#t_critical <- qt(0.975, df = n - 1) : better & reusable way to find critical t value for 95% confidence interval
qt(p = 0.025, df = 131, lower.tail = FALSE)
lowerCL <- ybar-(SEMy*1.978239) #This calculates:[lower limit = ybar - t*SE]
lowerCL
upperCL <- ybar+(SEMy*1.978239) #This calculates:[upper limit = ybar + t*SE]
upperCL

#complete reusable code for calculating confidence limits:
#{values <- data$Circumference..cm.
#values <- values[!is.na(values)]
#n <- length(values)
#ybar <- mean(values)
#s <- sd(values)
#SE <- s / sqrt(n)
#t_critical <- qt(0.975, df = n - 1)
#lowerCL <- ybar - t_critical * SE
#upperCL <- ybar + t_critical * SE
#lowerCL
#upperCL}

#make a histogram w/ the confidence limits
#hist(
#  orangedata$Circumference..cm., 
#  breaks = 30, 
#  main = "", 
#  xlab= "Circumference", 
#  ylim = c(0,25)
#  ) 
#abline(v=lowerCL, col="blue") #adds vertical line to histogram at lower confidence limit
#abline(v=upperCL, col="blue") #adds vertical line to histogram at upper confidence limit

#Compare quantiles of the data to the quantiles of a normal distribution
#qqnorm(orangedata$Circumference..cm.)
#qqline(orangedata$Circumference..cm., col = "red", lwd = 2)

# Select circumference values
circumference <- orangedata$Circumference..cm.

# Remove missing values
circumference <- circumference[!is.na(circumference)]

# Create normal Q-Q plot
qqnorm(
  circumference,
  main = "Normal Q-Q Plot of Tree Circumference",
  xlab = "Theoretical Normal Quantiles",
  ylab = "Observed Circumference"
)

# Add reference line
qqline(
  circumference,
  col = "red",
  lwd = 2
)

qqPlot(
  circumference,
  distribution = "norm", #Compares the quantiles of the data to the quantiles of a normal distribution
  envelope = 0.95, #Adds a 95% confidence envelope to the plot
  id = FALSE, #Prevents R from labeling outliers on the plot
  main = "Normal Q-Q Plot of Tree Circumference" #Adds a title to the plot
)