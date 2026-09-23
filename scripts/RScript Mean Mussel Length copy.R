#Imports Spreadsheet
data <- read.csv("Data/Byssal threads(sheet1).csv")

#Opens Spreadsheet In New Tab
View(data)

#Calculates Mean Length for Pristine Site
prmeanlen <- data %>%
  filter(Site=="Pristine") %>%
  summarise(mean(Length))

#Calculates Length Standard Deviation for Pristine Site
prstdlen <- data %>%
  filter(Site=="Pristine") %>%
  summarise(sd(Length))

#Calculates Mean Length for Polluted Site 
pomeanlen <- data %>%
  filter(Site=="Polluted") %>%
  summarise(mean(Length)) 

#Calculates Length Standard Deviation for Polluted Site
postdlen <- data %>%
  filter(Site=="Polluted") %>%
  summarise(sd(Length))
#Summary Function Above Generated Summary of Respective Data in a "mini-Table"

#Converts "mini-Table" Summary into Numeric Value
prmeanlen <- as.numeric(prmeanlen)
pomeanlen <- as.numeric(pomeanlen)

#Combines two values into single vector
y <- c(prmeanlen, pomeanlen)

#Creates site labels under bar graph - by making them into single line matrix 
names <- matrix(c("Pristine", "Polluted"))

#Creates data table that'll be used for the graph
graphdata <- data.frame(Site = names, Length = y)

#Creates the actual bar graph 
barplot(height = graphdata$Length, names.arg = graphdata$Site, 
        
#sets y-axis from 0 to 60        
        ylim = c(0,60), xlab = "Site", ylab = "Mean Length (cm)")

#Turns Pristine & Polluted SD from min-table into a number
prstdlen <- as.numeric(prstdlen)
postdlen <- as.numeric(postdlen)

#Combines the two numbers into single vector
stdev <- c(prstdlen, postdlen)

#Redraws graph & saves bar position into barploter
barploter <- barplot(height = graphdata$Length, names.arg = graphdata$Site, 
                     ylim = c(0,100), xlab = "Site", ylab = "Mean Length (cm)")

#Draws Error bars 
arrows(x0 = barploter, y0 = y - stdev,
       x1 = barploter, y1 = y + stdev,
       angle = 90, code = 3, length = 0.05)