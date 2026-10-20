#H_0: Polarized Wings of female Heliconius butterflies had no efect on mate attraction 
#H_A: Polarized Wings of female Heliconius butterflies had an effect on mate attraction 
#Test of choice: T-test 
read.csv("Data/Butterfly wings.csv")
data <- read.csv("Data/Butterfly wings.csv")

#subtract Difference Depolarized - Polarized 
Differences <- (data$Depolarized - data$Polarized)
Differences
data$Male

#Add collumn for Differences in Butterfly wings.csv /how to add collumn in csv file?data$Differences <- Differences
data$Differences <- (data$Depolarized - data$Polarized)
data

t.test(data$Polarized, data$Depolarized, paired = TRUE)

#Based on the t-test results, we can reject the null hypothesis and conclude that there is a significant difference in mate attraction between polarized and depolarized wings.ox
#Box
d <- data$Polarized - data$Depolarized
t.test(d, mu = 0)