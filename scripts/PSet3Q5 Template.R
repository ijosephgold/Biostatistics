# ============================================================
# 1. IMPORT THE DATA
# ============================================================

# Read the CSV file into R.
# Replace "Protein data.csv" with the name of another CSV file
# when you reuse this template.
protein_data <- read.csv("Data/Protein data.csv")


# ============================================================
# 2. INSPECT THE DATA
# ============================================================

# Display the first six rows.
# This helps you confirm that R imported the data correctly.
head(protein_data)

# Display the names of all columns.
# This tells you the exact names you must type in your code.
names(protein_data)

# Display the structure of the data.
# This shows whether each column is numeric, character, etc.
str(protein_data)


# ============================================================
# 3. CALCULATE THE RESPONSE VARIABLE
# ============================================================

# Calculate protein concentration for each mollusk.
#
# Protein is the numerator: the amount of protein measured.
# Weight is the denominator: the dry tissue mass.
#
# The $ symbol means "the column named..." inside a data frame.
# The <- symbol assigns/saves the result.
# The / symbol performs division.
protein_data$Concentration <-
  protein_data$Protein / protein_data$Weight

# Display the updated data frame so you can see the new column.
protein_data


# ============================================================
# 4. CREATE SHORTER OBJECT NAMES
# ============================================================

# Store the concentration column in an object called x.
# This is optional, but it makes later code shorter and easier to read.
x <- protein_data$Concentration

# Store the expected/reference value.
# For this assignment, the question gives us 0.50.
expected_mean <- 0.50


# ============================================================
# 5. CHECK THE DATA FOR APPROXIMATE NORMALITY
# ============================================================

# Create a histogram of the calculated concentrations.
#
# A histogram shows the distribution of the observations.
# We use it to look for strong skewness, unusual gaps, or outliers.
hist(
  x,
  main = "Histogram of Protein Concentration",
  xlab = "Protein concentration (mg protein/mg tissue)",
  col = "lightblue",
  border = "white"
)

# Create a normal Q-Q plot.
# If the points generally follow the line, the data are
# reasonably consistent with a normal distribution.
qqnorm(x)

# Add a reference line to the Q-Q plot.
qqline(x, col = "red", lwd = 2)

# Optional formal normality test.
# H0: the data come from a normal distribution.
# HA: the data do not come from a normal distribution.
shapiro.test(x)


# Histogram with more clearly defined bins
hist(
  x,
  breaks = 6,
  main = "Histogram of Protein Concentration",
  xlab = "Protein concentration",
  col = "lightblue",
  border = "white"
)

# Q-Q plot
qqnorm(x)
qqline(x, col = "red", lwd = 2)

# Boxplot
boxplot(
  x,
  main = "Boxplot of Protein Concentration",
  ylab = "Protein concentration",
  col = "lightblue"
)