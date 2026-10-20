#H_0: Polarized Wings of female Heliconius butterflies had no efect on mate attraction 
#H_A: Polarized Wings of female Heliconius butterflies had an effect on mate attraction 
#Test of choice: T-test 
read.csv("Data/Butterfly wings.csv")
data <- read.csv("Data/Butterfly wings.csv")

#subtract Difference Depolarized - Polarized 
#Add collumn for Differences
data$Differences <- (data$Depolarized - data$Polarized)
data

hist(data$Differences, main = "Histogram of Differences in Mate Attraction", xlab = "Differences (Depolarized - Polarized)", col = "lightblue", border = "black")
qqPlot(data$Differences, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Q-Q plot of Differences in Mate Attraction", col.lines = "red")
skewness(data$Differences)  
(data$Differences)
kurtosis(data$Differences)
t.test(data$Polarized, data$Depolarized, paired = TRUE)

#Based on the t-test results, we can reject the null hypothesis and conclude that there is a significant difference in mate attraction between polarized and depolarized wings.ox
#Box
d <- data$Polarized - data$Depolarized
t.test(d, mu = 0)