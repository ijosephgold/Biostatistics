#Imports Spreadsheet
data <- read.csv("Data/Byssal threads(sheet1).csv")

#Opens Spreadsheet In New Tab
View(data)

#Calculates Mean Thread.count for Pristine Site
prmeanthr <- data %>%
  filter(Site=="Pristine") %>%
  summarise(mean(Thread.count))

#Calculates Thread.count Standard Deviation for Pristine Site
prstdthr <- data %>%
  filter(Site=="Pristine") %>%
  summarise(sd(Thread.count))

#Calculates Mean Thread.count for Polluted Site 
pomeanthr <- data %>%
  filter(Site=="Polluted") %>%
  summarise(mean(Thread.count)) 

#Calculates Thread.count Standard Deviation for Polluted Site
postdthr <- data %>%
  filter(Site=="Polluted") %>%
  summarise(sd(Thread.count))
#Summary Function Above Generated Summary of Respective Data in a "mini-Table"

#Converts "mini-Table" Summary into Numeric Value
prmeanthr <- as.numeric(prmeanthr)
pomeanthr <- as.numeric(pomeanthr)

#Combines two values into single vector
y <- c(prmeanthr, pomeanthr)

#Creates site labels under bar graph - by making them into single line matrix 
names <- matrix(c("Pristine", "Polluted"))

#Creates data table that'll be used for the graph
graphdata <- data.frame(Site = names, Thread.count = y)

#Creates the actual bar graph 
barplot(height = graphdata$Thread.count, names.arg = graphdata$Site, 
        
        #sets y-axis from 0 to 60        
        ylim = c(0,60), xlab = "Site", ylab = "Mean Thread.count")

#Turns Pristine & Polluted SD from mini-table into a number
prstdthr <- as.numeric(prstdthr)
postdthr <- as.numeric(postdthr)

#Combines the two numbers into single vector
stdev <- c(prstdthr, postdthr)

#Redraws graph & saves bar position into barploter
barploter <- barplot(height = graphdata$Thread.count, names.arg = graphdata$Site, 
                     ylim = c(0,60), xlab = "Site", ylab = "Mean Thread.count")

#Draws Error bars 
arrows(x0 = barploter, y0 = y - stdev,
       x1 = barploter, y1 = y + stdev,
       angle = 90, code = 3, length = 0.05)