#Write data
nectar <- c(
  0.09, 0.13, 0.19, 0.25, 0.31, 0.35, 0.45, 0.56, 0.62,
  0.63, 0.75, 0.78, 1.16, 1.25, 1.34, 1.46, 2.30, 2.41,
  2.68, 2.80, 2.84, 3.28, 3.59, 4.27, 4.63, 5.40, 6.225,
  7.32, 9.30, 10.47, 10.67, 11.72, 13.59, 14.60, 24.63, 40.75
)

length(nectar)

#Examine
summary(nectar)

#create histogram 
par(mfrow = c(1, 2))

hist(
  nectar,
  main = "Original Nectar Volumes",
  xlab = "Nectar volume (µL)",
  col = "lightblue",
  border = "black"
)

boxplot(
  nectar,
  main = "Original Nectar Volumes",
  ylab = "Nectar volume (µL)",
  col = "lightblue"
)

par(mfrow = c(1, 1))
#data is strongly skewed to the right 

#try applying log transformation 
nectar_log <- log10(nectar)

summary(nectar_log)

par(mfrow = c(1, 2))

hist(
  nectar_log,
  main = "Log-Transformed Nectar Volumes",
  xlab = "log10(Nectar volume)",
  col = "lightgreen",
  border = "black"
)

boxplot(
  nectar_log,
  main = "Log-Transformed Nectar Volumes",
  ylab = "log10(Nectar volume)",
  col = "lightgreen"
)

par(mfrow = c(1, 1))


#comparing the two
par(mfrow = c(2, 2))

hist(
  nectar,
  main = "Original Histogram",
  xlab = "Nectar volume (µL)",
  col = "lightblue"
)

boxplot(
  nectar,
  main = "Original Boxplot",
  ylab = "Nectar volume (µL)",
  col = "lightblue"
)

hist(
  nectar_log,
  main = "Log-Transformed Histogram",
  xlab = "log10(Nectar volume)",
  col = "lightgreen"
)

boxplot(
  nectar_log,
  main = "Log-Transformed Boxplot",
  ylab = "log10(Nectar volume)",
  col = "lightgreen"
)

par(mfrow = c(1, 1))



#Reusable template for similar problems
# Enter the observations
data_name <- c(
# Put values here, separated by commas
)

# Check sample size
length(data_name)

# Examine original data
summary(data_name)

# Plot original data
#par sets or queries graphical parameters
par(mfrow = c(1, 2))

hist(
  data_name,
  main = "Original Data",
  xlab = "Measurement",
  col = "lightblue"
)

boxplot(
  data_name,
  main = "Original Data",
  ylab = "Measurement",
  col = "lightblue"
)

par(mfrow = c(1, 1))

# Apply log10 transformation
data_name_log <- log10(data_name)

# Examine transformed data
summary(data_name_log)

# Plot transformed data
par(mfrow = c(1, 2))

hist(
  data_name_log,
  main = "Log-Transformed Data",
  xlab = "log10(Measurement)",
  col = "lightgreen"
)

boxplot(
  data_name_log,
  main = "Log-Transformed Data",
  ylab = "log10(Measurement)",
  col = "lightgreen"
)

par(mfrow = c(1, 1))