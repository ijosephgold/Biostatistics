pen_data <- read.csv("data/African Penguin - Sheet1 (1).csv")

# 1. Check Normality & Homogeneity of Variance (Raw vs Transformed)
mod <- lm(X5 ~ Male.weight..kg., data = pen_data)
shapiro.test(residuals(mod))                  # Normality check on model residuals
bartlett.test(X5 ~ Male.weight..kg., data = pen_data)             # Bartlett's test for equal variance (Professor's primary test; p >= 0.05 = equal var)
fligner.test(X5 ~ Male.weight..kg., data = pen_data)              # Fligner-Killeen test (robust backup if non-normal)
library(moments); skewness(residuals(mod)); kurtosis(residuals(mod)) # Skew near 0 = symmetric; Kurt near 3 = normal tails
biostat_transform(pen_data)                         # Check if log10(Y) or sqrt(Y) fixes assumptions!

# 2. ASSUMPTION GRAPHS: Normal Q-Q Plot & Residuals vs. Fitted Values
#    WHY: Verifies normality (dots follow diagonal line) and equal spread across both groups.
library(car)
par(mfrow = c(1, 3))
res <- residuals(mod); hist(res, freq = FALSE, col = 'lightblue', border = 'white', main = 'Histogram')
qqPlot(res, col.lines = 'firebrick', main = 'QQ Plot w/ 95% Envelope') # Dots inside shaded gray/blue band = Normal!
plot(res ~ fitted(mod), pch = 19, xlab = 'Fitted Values', ylab = 'Residuals', main = 'Residuals vs Fitted'); abline(h = 0, lty = 2, col = 'red')
par(mfrow = c(1, 1))

# 3. Run Two-Sample Test Based on Assumptions
t.test(X5 ~ Male.weight..kg., data = pen_data, var.equal = FALSE) # Welch's t-test (unequal variances / safe default)
t.test(X5 ~ Male.weight..kg., data = pen_data, var.equal = TRUE)  # Student's t-test (equal variances)
wilcox.test(X5 ~ Male.weight..kg., data = pen_data)               # Mann-Whitney U test (if transformations fail)

# 4. PRIMARY FINAL GRAPH: Boxplot + Horizontal Jitter + Mean Diamond (Exact Wizard Plot!)
#    WHY THIS GRAPH: Shows full raw distribution, median/IQR, skewness, AND group mean all at once.
#    NOTE ON JITTER: `width = 0.15` nudges dots left/right ONLY so identical values don't overlap;
#                    vertical y-position is exact! `outlier.shape = NA` avoids drawing outliers twice.
library(ggplot2)
ggplot(pen_data, aes(x = factor(Male.weight..kg.), y = X5, fill = factor(Male.weight..kg.))) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +          # Horizontal jitter only
  stat_summary(fun = mean, geom = 'point', shape = 23, size = 3.5, fill = 'white') + # White Diamond = Mean
  labs(title = 'X5 by Male.weight..kg.',
       subtitle = 'Boxes = Median/IQR | White Diamond = Group Mean | Points = Raw Data',
       x = 'Male.weight..kg.', y = 'X5') +
  theme_minimal(base_size = 13) + theme(legend.position = 'none')

# 5. ALTERNATIVE GRAPH: Barplot of Group Means +/- SEM (Classic Exam / Publication Format)
ggplot(pen_data, aes(x = factor(Male.weight..kg.), y = X5, fill = factor(Male.weight..kg.))) +
  stat_summary(fun = mean, geom = 'bar', width = 0.55, color = 'black', alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', width = 0.18, linewidth = 0.8) +
  labs(title = 'X5 by Male.weight..kg. (Mean +/- SEM)', x = 'Male.weight..kg.', y = 'X5 (Mean +/- SEM)') +
  theme_classic(base_size = 13) + theme(legend.position = 'none')
