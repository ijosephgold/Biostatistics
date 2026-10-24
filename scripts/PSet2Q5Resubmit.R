#H_0: The limb lengths are the same for the hind leg and foreleg of deer
#H_A: The forelimbs are shorter than the hind limbs of deer
#Load the dataset into an R data frame
deer_data <- read.csv("Data/Deer legs.csv")
#View the data in a new tab with View()
View(deer_data)

#Perform the two-tailed paired t-test
#alternative="Greater" argument would be appropriate, but Profs hint said explicitly here to "use a two-tailed test" 
t.test(deer_data$Hindleg, deer_data$Foreleg, paired = TRUE, alternative = "two.sided")

#Conclusion: reject null hopothesis, hidlimbs and forelimbs differ significantly 
#Hidlegs are on average 3.3 longer than forelegs