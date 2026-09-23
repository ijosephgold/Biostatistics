data <- read.csv("Data/Orange trees FA26 - Sheet1.csv")
#View the data in a new tab with View()
View(data)
#Or just run it in the console
data

#create a binomial probability distribution of n = 132 trials
dist <- dbinom(x = 0:132, size = 132, p = 0.5)
dist
#Argument alternative =c("less")

#plot the distribution
barplot(dist, ylim= c(0,0.132), names.arg = 0:132, xlab = 'x', ylab = 'p')
#add color to the values of the distribution 0:9 
barplot(dist, ylim= c(0,0.132), names.arg = 0:132, xlab = 'x', ylab = 'p',
        col = c(rep("blue", 9), rep("gray", 11)))

p <- 2*dbinom(x=0:1, size = 132, p = 0.5)
sum(p) 


binom.test(1, 132, 0.5, alternative = c("less"))
