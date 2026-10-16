data <- read.csv("Data/Orange trees FA26 - Sheet1.csv")
#View the data in a new tab with View()
View(data)

total_trees <- nrow(data)

trees_with_fruit <- sum(data$Fruit.present == 1)
# view summarized relevant data 
print(paste("Total trees:", total_trees))
print(paste("Trees with fruit:", trees_with_fruit))

# now we can run one-tailed binomial test!
# x = the number of successes (132)
# n = the total number of trials (132)
# p = the hypothesized probability (0.95)
# alternative = "less" forces a one-tailed test checking if the true rate is lower than 95%
binom.test(x = trees_with_fruit, n = total_trees, p = 0.95, alternative = "less")

#Create a sequence of numbers, graph's scale?
x_values <- 0:total_trees

#Calculate theoretical binomial probabilities for each possible outcome in sequence, assuming nullH (p = 0.95) true
probabilities <- dbinom(x_values, size = total_trees, prob = 0.95)

#Draw bar-like plot (type="h" means histogram-like vertical lines)
plot(x_values, probabilities, type = "h", lwd = 5, col = "lightblue",
     main = "Binomial Distribution of Citrus Trees (n=132, p=0.95)",
     xlab = "Number of Trees with Fruit",
     ylab = "Probability")

#Add dashed red line showing where ACTUAL observed data landed (132)
abline(v = trees_with_fruit, col = "red", lwd = 2, lty = 2)

#Throw legend in the corner so graph is easy to read
legend("topleft", legend=c("Expected Probabilities", "Observed (132 trees)"), 
       col=c("lightblue", "red"), lwd=c(5, 2), lty=c(1, 2))
