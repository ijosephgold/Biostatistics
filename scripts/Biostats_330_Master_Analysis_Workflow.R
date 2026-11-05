# Biostats 330 — Master Analysis Workflow
#
# Purpose:
#   This is a decision guide and copy/paste template built from the R scripts,
#   datasets, transformations, assumption checks, and professor-style methods
#   already used in the Biostatistics repository.
#
# How to use this file:
#   1. Start at Section 0 every time.
#   2. Answer the simple QUESTIONS in the NOTE comments.
#   3. Run only the branch that matches your question.
#   4. Do not choose a test because it is familiar; choose it because the
#      question and study design require it.
#
# Important:
#   This file is a guide, not one giant script to run from top to bottom.
#   Replace words in ALL_CAPS before running a branch.


# ============================================================================
# 0. START EVERY ANALYSIS HERE
# ============================================================================

# NOTE — What is the first thing to do?
# Answer: Start with a clean RStudio session and make sure the project is the
# Biostatistics project. Then load packages. Installing and loading are not
# the same thing: install a package once; load it during each R session.

# Optional: use RStudio's Session > Restart R instead of deleting everything.
# rm(list = ls())

# The project files are organized as:
#   data/       CSV and Excel data
#   scripts/    question-specific scripts
#   output/     saved plots and results
#   biostats_toolkit.R  reusable functions
#   run_sync.R           GitHub synchronization; keep it separate

if (file.exists("biostats_toolkit.R")) {
  source("biostats_toolkit.R")
}

# This fallback loader makes the workflow usable even if the toolkit was not
# sourced. It installs only packages that are actually missing.
if (!exists("load_biostats_packages")) {
  load_biostats_packages <- function(install_missing = TRUE) {
    packages <- c(
      "tidyverse", "car", "moments", "DescTools", "ggplot2",
      "rstatix", "afex", "emmeans", "ggpubr", "agricolae",
      "lme4", "vegan", "readxl"
    )

    missing <- packages[!vapply(packages, requireNamespace,
                                 logical(1), quietly = TRUE)]

    if (length(missing) > 0 && install_missing) {
      install.packages(missing, repos = "https://cloud.r-project.org")
    }

    still_missing <- packages[!vapply(packages, requireNamespace,
                                      logical(1), quietly = TRUE)]
    if (length(still_missing) > 0) {
      stop("These packages are still missing: ",
           paste(still_missing, collapse = ", "))
    }

    invisible(lapply(packages, library, character.only = TRUE))
    message("Biostats packages are loaded.")
  }
}

load_biostats_packages(install_missing = TRUE)


# ============================================================================
# 1. BRING IN THE DATA AND PROVE THAT IT LOADED CORRECTLY
# ============================================================================

# NOTE — How should I import the data?
# Easiest answer: use file.choose() when you are unsure about the path.
# For a reproducible submission, replace it with a project-relative path such
# as "data/Orange trees FA26 - Sheet1.csv".

data <- read.csv(file.choose(), check.names = FALSE)

# For Excel files, use this instead:
# data <- readxl::read_excel(file.choose())

# NOTE — What should I look at immediately after importing?
# Ask these simple questions in order:
#   1. Did R create the object I expected?       -> exists("data")
#   2. Are the rows observations and columns variables? -> head(data)
#   3. Which columns are numeric, categorical, or IDs? -> str(data)
#   4. Are there missing values?                 -> colSums(is.na(data))
#   5. Are numbers accidentally characters/factors? -> str(data)
#   6. Do the minimum and maximum values make biological sense? -> summary(data)

exists("data")
dim(data)
names(data)
head(data)
tail(data)
str(data)
summary(data)
colSums(is.na(data))

# A compact inspection helper used throughout your scripts:
inspect_data <- function(x) {
  cat("Rows:", nrow(x), " Columns:", ncol(x), "\n\n")
  print(names(x))
  cat("\nStructure:\n")
  print(str(x))
  cat("\nSummary:\n")
  print(summary(x))
  cat("\nMissing values by column:\n")
  print(colSums(is.na(x)))
  invisible(x)
}

inspect_data(data)

# NOTE — What if the file will not load?
#   - If the error says the file cannot be found, use file.choose().
#   - If the columns are in one giant column, check the delimiter and use
#     read.csv2() or read.delim() only if the file actually uses that format.
#   - If a column name contains spaces, use data[["Column name"]] or rename it.
#   - In this repository, use lowercase data/ rather than Data/ when possible.


# ============================================================================
# 2. WRITE THE QUESTION BEFORE RUNNING A TEST
# ============================================================================

# NOTE — What am I trying to compare or relate?
# Choose the first statement that is true:
#
#   A. One numeric measurement compared with a fixed value?
#      -> one-sample t-test or its nonparametric alternative.
#
#   B. One proportion/count compared with a stated probability or ratio?
#      -> exact binomial, Pearson chi-square, or G-test.
#
#   C. Two categorical variables, such as color by location?
#      -> contingency table; chi-square/G-test or Fisher's exact test.
#
#   D. Two numeric measurements from two independent groups?
#      -> independent two-sample t-test, Welch t-test, or Wilcoxon.
#
#   E. Two measurements from the same individual, or matched pairs?
#      -> paired t-test, paired Wilcoxon, or paired sign test.
#
#   F. A numeric measurement across three or more independent groups?
#      -> one-way ANOVA, Welch ANOVA, or Kruskal-Wallis.
#
#   G. The same subjects measured under three or more conditions?
#      -> repeated-measures ANOVA, with subject ID and within-subject factor.
#
#   H. The relationship between two numeric variables?
#      -> correlation; use regression too when a prediction equation is wanted.

# Write these before the code:
#   H_0: [what no difference/no relationship would mean]
#   H_A: [what difference/relationship you are looking for]
#   alpha <- 0.05

alpha <- 0.05

# NOTE — How do I choose "less", "greater", or "two.sided"?
# Ask: Did the question predict a direction before seeing the data?
#   - Yes, lower -> alternative = "less"
#   - Yes, higher -> alternative = "greater"
#   - No specific direction, or the professor says two-tailed -> "two.sided"
# Do not change the direction after seeing the result.


# ============================================================================
# 3. IDENTIFY THE DESIGN — THE DECISION THAT CONTROLS EVERYTHING AFTER IT
# ============================================================================

# NOTE — Are the observations independent?
# Ask: Could one observation's measurement influence another observation's
# measurement, or were both measurements taken from the same subject/tree/deer?
#   - Same individual or deliberately matched units -> paired/repeated design.
#   - Different subjects/trees/plots with no matching -> independent design.
#   - If unsure, look at the study description, not at the p-value.
# Independence is determined by how the data were collected; R cannot repair
# a design that is not independent.

# NOTE — Is the response numeric or categorical?
#   - Numeric: length, mass, height, temperature, concentration, diversity.
#   - Categorical/count: red/black, rural/industrial, observed category counts.
#   - ID columns such as Individual, Subject, Tree number, or Quadrat are not
#     outcomes. Do not test an ID column.


# ============================================================================
# 4. DESCRIPTIVE GRAPHS BEFORE FORMAL ASSUMPTION TESTS
# ============================================================================

# NOTE — What should I graph?
#   - One numeric variable: histogram, boxplot, Q-Q plot.
#   - Numeric response by group: group boxplot and/or jittered points.
#   - Two numeric variables: scatterplot with a fitted line.
#   - Counts/categories: barplot or contingency-table barplot.
# The graph answers the practical question: does the data look roughly
# symmetric, does one group spread much more than another, and are there
# outliers or impossible values?

# Replace RESPONSE, GROUP, X, and Y with real column names.
# hist(data$RESPONSE, main = "", xlab = "RESPONSE")
# boxplot(RESPONSE ~ GROUP, data = data, xlab = "GROUP", ylab = "RESPONSE")
# plot(data$X, data$Y, pch = 19, xlab = "X", ylab = "Y")
# abline(lm(Y ~ X, data = data), col = "blue", lwd = 2)


# ============================================================================
# 5. ASSUMPTIONS — SIMPLE QUESTIONS AND WHAT EACH ANSWER MEANS
# ============================================================================

normality_report <- function(x, label = "values", make_plots = TRUE) {
  x <- x[is.finite(x)]
  cat("\n--- Normality report:", label, "---\n")
  cat("n =", length(x), "\n")

  if (length(x) >= 3 && length(x) <= 5000) {
    print(shapiro.test(x))
  } else {
    cat("Shapiro-Wilk skipped: n must be between 3 and 5000.\n")
  }

  if (length(x) >= 3) {
    cat("Skewness:", moments::skewness(x), "\n")
    cat("Kurtosis:", moments::kurtosis(x), "\n")
  }

  if (make_plots) {
    hist(x, main = paste("Histogram:", label), xlab = label)
    car::qqPlot(x, distribution = "norm",
                main = paste("Q-Q plot:", label))
  }
  invisible(x)
}

# NOTE — What does Shapiro-Wilk test?
# H_0: the values are consistent with a normal distribution.
# If p <= 0.05, there is evidence against normality.
# If p > 0.05, you did not detect evidence against normality; this does not
# prove the data are perfectly normal. Always inspect the histogram and Q-Q
# plot too, especially with small samples.

# NOTE — Which values should I test for normality?
#   - One-sample test: the response values.
#   - Independent two-group test: the response within each group.
#   - Paired test: the paired differences, not the two raw columns separately.
#   - ANOVA: the model residuals, not only the raw response.
#   - Regression: the model residuals, plus the residual-vs-fitted plot.
#   - Repeated measures ANOVA: the model residuals and the study design.

# Example:
# normality_report(data$RESPONSE, "RESPONSE")
# normality_report(data$RESPONSE[data$GROUP == "Group 1"], "Group 1")

# NOTE — How do I check equal variances?
#   1. Look at side-by-side boxplots first.
#   2. If the distributions are reasonably normal, use Bartlett as in your
#      professor's scripts.
#   3. If normality is doubtful, use Levene's test as the more robust check.
# H_0 for both: the group variances are equal.
#   - p > 0.05: no evidence of unequal variances detected.
#   - p <= 0.05: treat the variances as unequal.

# Bartlett example:
# bartlett.test(RESPONSE ~ GROUP, data = data)
# Levene example:
# car::leveneTest(RESPONSE ~ GROUP, data = data)


# ============================================================================
# 6. TRANSFORM ONLY WHEN THE CHECKS SHOW A REASON
# ============================================================================

# NOTE — Do I need a transformation?
# Ask: Are skewness, the Q-Q plot, or the residual plot bad enough that the
# planned parametric analysis is not reasonable?
#   - No -> keep the original response.
#   - Yes -> try the transformation that matches the data, then repeat the
#     same graph and assumption checks.
#
# Transformations used in your scripts:
#   - sqrt(RESPONSE): useful for nonnegative counts and moderate right skew;
#     works when zeros are present.
#   - log(RESPONSE): useful for positive right-skewed measurements; it cannot
#     be used on zero or negative values.
#   - Do not transform just to make a p-value cross 0.05.
#
# Example:
# data$RESPONSE_sqrt <- sqrt(data$RESPONSE)
# data$RESPONSE_log <- log(data$RESPONSE)  # only if every value is > 0
# normality_report(data$RESPONSE_sqrt, "sqrt(RESPONSE)")
# normality_report(data$RESPONSE_log, "log(RESPONSE)")

# NOTE — Which version should I keep?
# Keep the version with the most reasonable biological scale and the best
# overall diagnostics: Q-Q plot, histogram, residual pattern, spread, and
# skewness. Shapiro-Wilk alone does not choose the transformation.
# The final statistical test must use the chosen version consistently.


# ============================================================================
# 7. BRANCH A — ONE NUMERIC SAMPLE AGAINST A FIXED VALUE
# ============================================================================

# Used in your work for questions such as Buckwheat density compared with 4.

# NOTE — Is the question "is the population mean different from VALUE"?
#   - Yes -> one-sample t-test if the response is reasonably normal.
#   - If the response remains strongly nonnormal, ask whether the professor
#     wants a one-sample sign/Wilcoxon method.

# H_0: mean(RESPONSE) = VALUE
# H_A: mean(RESPONSE) != VALUE  # change direction only if predicted

# t.test(data$RESPONSE, mu = VALUE, alternative = "two.sided")

# Report: sample mean, confidence interval, t, df, p-value, and the decision
# about H_0. If significant, state the direction using the sample mean.


# ============================================================================
# 8. BRANCH B — ONE PROPORTION OR ONE CATEGORICAL RATIO
# ============================================================================

# B1. Exact binomial test — used for the Orange tree fruit question.

# NOTE — Is each observation a success/failure and is the null a probability?
#   - Yes -> use binom.test().
#   - Count successes with sum(), count trials with nrow() or the stated n.

# total <- nrow(data)
# successes <- sum(data$SUCCESS_COLUMN == 1, na.rm = TRUE)
# binom.test(successes, total, p = 0.95, alternative = "less")

# NOTE — What does the alternative mean here?
# "less" asks whether the true success rate is below the stated rate.
# "greater" asks whether it is above the stated rate.
# "two.sided" asks whether it differs in either direction.

# B2. Goodness-of-fit chi-square/G-test — used for bee sex ratio and the
# dihybrid 9:3:3:1 question.

# NOTE — Is one categorical variable being compared with a theoretical ratio?
#   - Yes -> create observed counts and expected counts in the same order.
#   - Expected counts = total count * expected proportions.
#   - df = number of categories - 1.

# observed <- matrix(c(COUNT_1, COUNT_2, COUNT_3), nrow = 1)
# expected <- matrix(sum(observed) * c(P_1, P_2, P_3), nrow = 1)

# Pearson chi-square, in the manual style used in your scripts:
# X2 <- sum((observed - expected)^2 / expected)
# p_X2 <- pchisq(X2, df = length(observed) - 1, lower.tail = FALSE)

# G-test, in the manual style used in your scripts:
# contribution <- ifelse(observed == 0, 0,
#                        observed * log(observed / expected))
# G2 <- 2 * sum(contribution)
# p_G <- pchisq(G2, df = length(observed) - 1, lower.tail = FALSE)

# NOTE — What do I conclude?
#   - p <= alpha: reject H_0; the observed counts do not fit the proposed
#     ratio closely enough.
#   - p > alpha: fail to reject H_0; there is not enough evidence that the
#     observed counts differ from the proposed ratio.
# Do not say "the null is proven true."

# B3. Intrinsic/Poisson goodness-of-fit — used for count frequencies such as
# oranges per tree and the Midterm Q3/Mites-and-midges style problem.

# NOTE — Is the question whether count frequencies follow a Poisson model?
#   - Yes -> estimate lambda from the observed counts, calculate expected
#     Poisson probabilities, convert them to expected counts, and combine rare
#     high-count categories (for example, a final "4+" category).
#   - When lambda is estimated from the same data, use df = final number of
#     categories - 2, as in the professor's notes.
#   - Then compare the observed and expected counts with the required G-test
#     or chi-square calculation.


# ============================================================================
# 9. BRANCH C — TWO CATEGORICAL VARIABLES: CONTINGENCY TABLES
# ============================================================================

# Used for Ladybird color by location, Foraging behavior by substrate, and
# bird-guild style data.

# NOTE — Does the question ask whether two categories are associated?
#   - Yes -> make a table with rows = one categorical variable and columns =
#     the other. Do not use raw individual IDs as categories.

# If the data are already summarized counts:
# observed <- matrix(c(COUNT_A, COUNT_B, COUNT_C, COUNT_D),
#                    nrow = 2, byrow = TRUE)
# rownames(observed) <- c("ROW_1", "ROW_2")
# colnames(observed) <- c("COL_1", "COL_2")

# If each row is an observation with a count column:
# observed <- xtabs(COUNT_COLUMN ~ ROW_COLUMN + COL_COLUMN, data = data)

# addmargins(observed)
# expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
# expected

# NOTE — Are the expected counts large enough for the planned chi-square/G-test?
#   - Yes -> use the professor-required Pearson chi-square or G-test.
#   - No, or are there structural zeros/sparse cells? -> use Fisher's exact
#     test for a small table, or the simulated Fisher procedure already used
#     in PSet3Q3.R. State the limitation.
#
# Important: an observed zero does not automatically mean "use nonparametric."
# G-tests are for categorical counts; the issue is sparse/invalid expected
# counts, not whether one observed cell happens to equal zero.

# Pearson chi-square:
# chisq.test(observed, correct = FALSE)

# G-test:
# DescTools::GTest(observed, correct = "none")

# Fisher/simulation branch used in PSet3Q3.R:
# fisher.test(observed, simulate.p.value = TRUE, B = 1000)

# Optional 2 x 2 association index used in the professor-style work:
# association_index <- ((observed[1, 1] + observed[2, 2]) -
#                       (observed[1, 2] + observed[2, 1])) / sum(observed)

# NOTE — What should the graph show?
#   - Counts: side-by-side bars for raw observed counts.
#   - Proportions: stacked/fill bars when comparing composition within rows.
# The graph describes the pattern; the test decides whether the association
# is stronger than expected by chance.


# ============================================================================
# 10. BRANCH D — TWO INDEPENDENT NUMERIC GROUPS
# ============================================================================

# Used in your work for examples such as Microalgae concentration and the
# Barnacle-style two-sample script.

# NOTE — Are the two groups independent?
#   - Yes -> continue here.
#   - No, same individuals/paired units -> go to Branch E.

# First inspect each group:
# boxplot(RESPONSE ~ GROUP, data = data)
# normality_report(data$RESPONSE[data$GROUP == "GROUP_1"], "GROUP_1")
# normality_report(data$RESPONSE[data$GROUP == "GROUP_2"], "GROUP_2")
# car::leveneTest(RESPONSE ~ GROUP, data = data)

# NOTE — Which test do I run?
#   1. Both groups look reasonably normal and spread is similar:
#      -> pooled two-sample t-test, var.equal = TRUE.
#   2. Groups look reasonably normal but spread is different:
#      -> Welch t-test, var.equal = FALSE.
#   3. Strong nonnormality remains after a sensible transformation or the
#      response is ordinal/outlier-dominated:
#      -> Wilcoxon rank-sum test, if permitted by the assignment.

# pooled t-test:
# t.test(RESPONSE ~ GROUP, data = data, var.equal = TRUE)

# Welch t-test:
# t.test(RESPONSE ~ GROUP, data = data, var.equal = FALSE)

# Wilcoxon rank-sum:
# wilcox.test(RESPONSE ~ GROUP, data = data)

# NOTE — Why does var.equal matter?
# var.equal = TRUE assumes the two population variances are equal.
# var.equal = FALSE is Welch's method and is safer when the variances differ.
# It does not mean the data are nonparametric; it is still a t-test for means.


# ============================================================================
# 11. BRANCH E — TWO PAIRED NUMERIC MEASUREMENTS
# ============================================================================

# Used for deer hindleg vs foreleg, butterfly wings, and paired cone/condition
# examples. The pairing must be preserved by row.

# NOTE — Is row 1's first measurement matched to row 1's second measurement?
#   - Yes -> paired analysis.
#   - No -> independent-groups analysis instead.

# differences <- data$MEASUREMENT_1 - data$MEASUREMENT_2
## normality_report(differences, "paired differences")

# NOTE — What is the assumption decision based on?
# The normality of the DIFFERENCES, not whether each raw column is normal.

# If differences are reasonably normal:
# t.test(data$MEASUREMENT_1, data$MEASUREMENT_2,
#        paired = TRUE, alternative = "two.sided")

# If differences remain strongly nonnormal:
# wilcox.test(data$MEASUREMENT_1, data$MEASUREMENT_2,
#             paired = TRUE, alternative = "two.sided")

# If the professor specifically asks for a paired sign test, use the sign of
# each nonzero difference and report that method explicitly.

# NOTE — What if the assignment gives an explicit direction?
# Follow the assignment's direction. For example, PSet2Q5 explicitly called
# for a two-tailed paired t-test even though the biological wording suggested
# a directional expectation.


# ============================================================================
# 12. BRANCH F — THREE OR MORE INDEPENDENT NUMERIC GROUPS
# ============================================================================

# Used for Erodium, Aphid length by gall, diversity by location, Sagebrush,
# and similar ANOVA questions.

# NOTE — Are there three or more independent groups?
#   - Yes -> one-way ANOVA family.
#   - No -> use the two-group branch instead.

# model <- lm(RESPONSE ~ GROUP, data = data)
# anova(model)
# summary(model)

# Check ANOVA assumptions on residuals:
# residuals_model <- resid(model)
# hist(residuals_model, main = "", xlab = "Residuals")
# car::qqPlot(residuals_model, distribution = "norm")
# shapiro.test(residuals_model)
# moments::skewness(residuals_model)
# moments::kurtosis(residuals_model)
# bartlett.test(residuals_model ~ data$GROUP)

# NOTE — What is the simple decision?
#   1. Residuals look acceptable and variances are similar:
#      -> ordinary ANOVA; if significant, use Tukey HSD.
#   2. Residuals are acceptable but variances differ:
#      -> Welch one-way ANOVA; use Games-Howell for pairwise follow-up if
#         permitted by the assignment.
#   3. Strong nonnormality remains after a sensible transformation:
#      -> Kruskal-Wallis, if this is the method allowed for the question.
#
# A significant ANOVA answers only "at least one group differs." It does not
# identify which groups differ until the post-hoc comparison is run.

# ordinary ANOVA and Tukey:
# anova(model)
# TukeyHSD(aov(model))

# Welch ANOVA and Games-Howell:
# oneway.test(RESPONSE ~ GROUP, data = data, var.equal = FALSE)
# rstatix::games_howell_test(data, RESPONSE ~ GROUP)

# Kruskal-Wallis:
# kruskal.test(RESPONSE ~ GROUP, data = data)

# NOTE — What if there is a random factor?
# Your ANOVA script includes Genet as a random factor. Do not replace that
# design with an ordinary one-way ANOVA without checking the study structure.
# The relevant model pattern is:
# lme4::lmer(RESPONSE ~ 1 + (1 | RANDOM_FACTOR), data = data)


# ============================================================================
# 13. BRANCH G — REPEATED-MEASURES ANOVA
# ============================================================================

# Used for Stress/BPM and Coral/FvFm questions.

# NOTE — Did the same subject/coral/tree receive multiple conditions or days?
#   - Yes -> repeated-measures design.
#   - No -> ordinary one-way ANOVA.

# The data must be in long format: one row per subject-condition measurement.
# Required columns are:
#   DV          numeric response
#   ID          subject/coral/tree identifier
#   WITHIN     repeated condition/day/test factor

# rm_model <- afex::aov_ez(
#   data = data,
#   dv = "DV",
#   id = "ID",
#   within = "WITHIN"
# )
# rm_model
# emmeans::emmeans(rm_model, ~ WITHIN)

# NOTE — What does the ID column do?
# It tells R which rows belong to the same individual. Without it, R cannot
# model the repeated measurements correctly.

# NOTE — What comes after a significant repeated-measures ANOVA?
# Use emmeans to compare the conditions and explain which condition means
# differ. Plot the means with uncertainty bars, but do not use the graph alone
# to decide significance.


# ============================================================================
# 14. BRANCH H — CORRELATION AND LINEAR REGRESSION
# ============================================================================

# Used for penguin beak width/depth, environmental variables, ants, transects,
# and stream biomass/species richness.

# NOTE — Is the question about association between two numeric variables?
#   - Yes -> make a scatterplot first.

# plot(data$X, data$Y, pch = 19, xlab = "X", ylab = "Y")

# NOTE — Does the scatterplot look approximately linear, with no extreme
# outlier controlling the pattern?
#   - Yes -> Pearson correlation is reasonable.
#   - No, but the pattern is monotonic -> Spearman correlation may be better.
#   - No clear pattern -> do not force a correlation just because both columns
#     are numeric.

# Pearson:
# cor.test(data$X, data$Y, method = "pearson")

# Spearman:
# cor.test(data$X, data$Y, method = "spearman")

# NOTE — Do I need a linear model too?
#   - If you only need strength/direction of association -> cor.test may be
#     enough.
#   - If you need a fitted equation, slope, prediction, or residual diagnostics
#     -> fit lm().

# model <- lm(Y ~ X, data = data)
# summary(model)
# plot(model, which = 1)  # residuals vs fitted: random cloud around 0
# plot(model, which = 2)  # normal Q-Q plot of residuals
# shapiro.test(resid(model))

# NOTE — What do I look for in the residual-vs-fitted plot?
#   - Random cloud around zero -> linearity/spread look reasonable.
#   - Curved pattern -> relationship may not be linear.
#   - Funnel shape -> spread changes with fitted values.
#   - One point far away -> possible influential outlier.

# For many correlations, as in PSet6Q2, adjust the p-values:
# variables <- c("X1", "X2", "X3", "X4")
# pairs <- combn(variables, 2, simplify = FALSE)
# results$p_adjusted <- p.adjust(results$p, method = "holm")


# ============================================================================
# 15. WIDE VERSUS LONG DATA — ONLY RESHAPE WHEN THE METHOD NEEDS IT
# ============================================================================

# NOTE — Is one row already one subject and separate columns are conditions?
#   - Paired t-test can often use the wide columns directly.
#   - Repeated-measures ANOVA and grouped plots usually need long format.

# Example wide -> long using the style already used in your scripts:
# long_data <- data %>%
#   tidyr::pivot_longer(
#     cols = c(CONDITION_1, CONDITION_2),
#     names_to = "Condition",
#     values_to = "Response"
#   )

# Older scripts use gather(); pivot_longer() is the newer equivalent:
# long_data <- data %>%
#   tidyr::gather(key = "Condition", value = "Response",
#                 CONDITION_1, CONDITION_2, na.rm = TRUE)

# NOTE — What should stay the same after reshaping?
# The actual measurements and subject IDs. Reshaping changes the layout, not
# the observations. Check nrow(), the IDs, and the first few rows afterward.


# ============================================================================
# 16. FINAL RESULTS CHECKLIST
# ============================================================================

# Before writing the answer, ask:
#   1. Did I use the test that matches the design?
#   2. Did I state H_0 and H_A in words?
#   3. Did I identify the correct test statistic, df, and p-value?
#   4. Did I use alpha = 0.05 unless the question says otherwise?
#   5. Did I say reject or fail to reject H_0, not "accept H_0"?
#   6. If significant, did I state the direction or which groups differ?
#   7. If ANOVA was significant, did I run the required post-hoc comparison?
#   8. If I transformed the response, did I say which scale the test used?
#   9. Did I report biologically meaningful means/medians, differences, or r?
#  10. Did I save the graph and keep the script reproducible?

# Result sentence templates:
#
# #notes — p-value decision:
# "Because p = [p] is [less/greater] than alpha = 0.05, I [reject/fail to
# reject] H_0."
#
# #notes — numeric comparison:
# "The [test] showed that [response] was [higher/lower/different] in [group]
# than [comparison], [statistic] = [value], df = [value], p = [value]."
#
# #notes — correlation:
# "There was a [positive/negative] [weak/moderate/strong] association between
# X and Y, r = [r], p = [p]. This describes association, not causation."
#
# #notes — non-significant result:
# "The result was not statistically significant, so I fail to reject H_0. This
# means the data did not provide enough evidence for the proposed difference;
# it does not prove that the groups are identical."


# ============================================================================
# 17. QUESTION-TO-DATA MAP FROM THE CURRENT BIOSTATISTICS REPOSITORY
# ============================================================================

# Categorical/count questions:
#   PSet2Q1Resubmit.R       -> bee sex ratio, manual X2 and G2
#   PSet2Q2Resubmit.R       -> dihybrid 9:3:3:1 G-test
#   PSet2Q4Resubmit.R       -> Orange tree fruit proportion, exact binomial
#   Midterm Q3.R            -> Poisson-style goodness-of-fit
#   PSet3Q2.R               -> Ladybird color by location, G-test
#   PSet3Q3.R               -> Foraging behavior by substrate, Fisher simulation
#   contingency tables2.R   -> contingency-table practice
#
# One-sample, paired, or two-group questions:
#   Midterm Q2.R            -> Buckwheat density versus a fixed value
#   PSet2Q5Resubmit.R       -> paired deer leg measurements
#   PSet3Q5.R               -> Protein-related t-test workflow
#   PSet4 Q3.R              -> Urchin measurements by habitat
#   PSet4 Q4.R              -> Microalgae two-group comparison and transformation
#   PSet4 Q5.R              -> paired Yew cone counts, Wilcoxon approach
#   PSet4Q2.R               -> paired butterfly wing measurements
#   Two sample t-test.R     -> independent groups, transformation, Wilcoxon
#
# ANOVA/repeated-measures questions:
#   ANOVA.R                -> Sagebrush, Genet/random-factor examples
#   PSet5Q2.R              -> Erodium treatment comparison
#   PSet5Q3.R              -> Aphid length by gall
#   PSet5Q5.R              -> Diversity by location, transformations/Kruskal
#   Repeated measures ANOVA.R -> Stress/BPM by repeated test
#   QSet5Q4.R / PSet5Q4 fancy.R -> Coral FvFm by repeated day/temperature
#
# Correlation/regression questions:
#   PSet6Q1.R              -> penguin beak width/depth
#   PSet6Q2.R              -> environmental-variable correlation matrix
#   Correlation.R          -> ant and transect correlations
#   Regression.R           -> stream biomass/species-richness regression
#
# Assumption/transformation/tool files:
#   Data transformations.R, OrangeTree_Circumference_Script.R,
#   assumption_report.R, assumption_report_interactive.R, inspect_df.R,
#   biostats_toolkit.R, export_plot.R

# Repository note:
#   Two existing scripts reference data files not currently present in the
#   current repository snapshot: Chaparral data.csv and Barnacle data.csv.
#   Several older scripts use Data/ with a capital D; use data/ in the current
#   repository or use file.choose() until those paths are standardized.

