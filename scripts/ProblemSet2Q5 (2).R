#Load the dataset into an R data frame
deer_data <- read.csv("Data/Deer legs.csv")
#View the data in a new tab with View()
View(deer_data)
#Or just run it in the console
data

#Perform the two-tailed paired t-test
# The alternative="two.sided" argument is R's default, but it is explicitly written here to match your prompt's hint.
t.test(deer_data$Hindleg, deer_data$Foreleg, paired = TRUE, alternative = "two.sided")
