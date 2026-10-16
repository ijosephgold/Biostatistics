#For a future problem, you mainly change the items in uppercase:
# ============================================================
# ONE-SAMPLE TEST TEMPLATE
# ============================================================

# Import the data
data <- read.csv("FILE_NAME.csv")

# Inspect the data
head(data)
names(data)
str(data)

# Calculate the response variable
# Replace NUMERATOR_COLUMN and DENOMINATOR_COLUMN
# with the actual column names in your file.
data$response <- data$NUMERATOR_COLUMN /
  data$DENOMINATOR_COLUMN

# Store the response variable
x <- data$response

# Enter the expected/reference value from the question
expected_mean <- REFERENCE_VALUE

# Examine normality
hist(
  x,
  main = "Histogram of Response Variable",
  xlab = "Response variable",
  col = "lightblue",
  border = "white"
)

qqnorm(x)
qqline(x, col = "red", lwd = 2)

# Optional normality test
shapiro.test(x)

# Run the one-sample test
test_result <- t.test(
  x,
  mu = expected_mean,
  alternative = "two.sided",
  conf.level = 0.95
)

# View the test result
test_result

# Make the boxplot
boxplot(
  x,
  main = "Boxplot of Response Variable",
  ylab = "Response variable",
  col = "lightblue"
)

# Add the expected/reference mean
abline(
  h = expected_mean,
  col = "red",
  lwd = 2,
  lty = 2
)

# Add the 95% confidence limits
abline(
  h = test_result$conf.int,
  col = "darkgreen",
  lwd = 2,
  lty = 3
)

# Add a legend
legend(
  "topright",
  legend = c(
    "Expected mean",
    "95% confidence limits"
  ),
  col = c("red", "darkgreen"),
  lty = c(2, 3),
  lwd = 2
)

# Make an automatic conclusion
if (test_result$p.value < 0.05) {
  print("Reject H0: the mean differs significantly from the expected value.")
} else {
  print("Fail to reject H0: no significant difference was detected.")
}