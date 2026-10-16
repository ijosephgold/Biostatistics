# Import CSV file
protein_data <- read.csv("../Data/Protein data.csv")
#view data
protein_data
# Inspect the data
head(protein_data)
names(protein_data)
str(protein_data)

#Calculate protein concentration
# The $ symbol means "the column named..." inside a data frame.
protein_data$Concentration <- 
  protein_data$Protein / protein_data$Weight

#viewing calculated concentrations
protein_data

#checking normality
hist( 
  protein_data$Concentration, 
  main = "Histogram of Protein Concentration", 
  xlab = "Protein Concentration (mg protein/mg tisue", 
  col = "lightblue", 
  border = "white" 
  )

qqnorm(protein_data$Concentration)
qqline(protein_data$Concentration, col = "red", lwd = 2)

#Shapiro Wilk test
shapiro.test(protein_data$Concentration)

#Execute Hypothesis test
expected <- 0.50

test_result <- t.test( 
  protein_data$Concentration, 
  mu = expected, 
  alternative = "two.sided", 
  con.level = 0.95 
  ) 

test_result

#Boxplot w expected mean & confidence limits
boxplot( 
  protein_data$Concentration, 
  main = "Protein Concentration in Marine Bivalves", 
  ylab = "Protein Concentration (mg protein/mg tissue)", 
  col = "lightblue" 
  ) 

#add expected mean
#albine adds a straight line to a graph
#h means horizontal 
abline( 
  h = expected, 
  col = "red", 
  lwd = 2, 
  lty = 2 
  ) 

#adding 95% confidence limits 
abline( 
  h = test_result$conf.int, 
  col = "darkgreen", 
  lwd = 2, 
  lty = 3 
  ) 

legend( 
  "topright", 
  legend = c("Expected mean = 0.50", "95% confidence limits"), 
  col = c("red", "darkgreen"), 
  lty = c(2, 3), 
  lwd = 2 
  )