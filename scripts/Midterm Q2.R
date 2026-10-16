#Since the Q asks if density "differs", use a 2 tailed test 
#H_0: Average Buckwheat density is 4 plants per 25m^2 
#H_A: Average Buckwheat density differs from 4 plants per 25m^2
buckwheat <- read.csv("C:/Biostatistics/data/Buckwheat data.csv")

names(buckwheat)
str(buckwheat)

buckwheat_density <- buckwheat$Density

#Ttest does sample mean differ from 4
#t.test(buckwheat_density, mu = 4, alternative = "two.sided")
test <- t.test(
  buckwheat_density,
  mu = 4,
  conf.level = 0.95
)

test

#Boxplot of Buckwheat density
boxplot(
  buckwheat_density,
  notch = TRUE,
  ylab = "Eriogonum plants per 25m^2",
  main = "Buckwheat Density",
  col = "lightgreen"
)
#h adds horizontal H_0 line to boxplot
abline(
  h = 4,
  col = "red",
  lty = 2,
  lwd = 2
)

#p statistic < than 0.05, reject H_0. Buckwheat density is significantly greater than 4/25m^2