# My R Scripts

## ANOVA

```r
#ANOVA

#load libraries
library(tidyverse)
library(car) #for QQ-plot
library(moments) #for skewness and kurtosis)
library(agricolae) #for Tukey's HSD
library(ggplot2)
library(lme4) #for random factors
library(vegan) #for perMANOVA

#Model I ANOVA
#test California Sagebrush heights in different soil types

#load data
sagedata <- read.csv("Data/Sagebrush data.csv")
sagedata
inspect_df(sagedata) #check for missing data

#fit data to a linear model with soil = x and height = y
sagedata$Soil <- as.factor(sagedata$Soil)
lm <- lm(Height~Soil, data=sagedata) #takes the data as lm(y~x)
lm
anova(lm) #provides the ANOVA table
summary(lm) #provides extra details about the model

#check the assumptions
#recall, we want to look at the residuals to assess normality
lmres <- resid(lm)
qqPlot(lmres, "norm")
hist(lmres, xlab = "Height residuals", main = "")
boxplot(lmres, ylab = "Height residuals")
shapiro.test(lmres) #P > 0.05
skewness(lmres) 
kurtosis(lmres)
#next test equality of variances
bartlett.test(Height~Soil, data=sagedata) #P > 0.05

#or make a boxplot
Sterile <- sagedata %>%
  filter(Soil == "Sterile")
SandDune <- sagedata %>%
  filter(Soil == "Sand Dune")
Mountains <- sagedata %>%
  filter(Soil == "Mountains")
Chaparral <- sagedata %>%
  filter(Soil =="Chaparral")
boxplot(Sterile$Height, SandDune$Height, Mountains$Height, Chaparral$Height,
        names = c("Sterile", "Sand Dune", "Mountains", "Chaparral"), 
        xlab = "Soil type", ylab = "Height (cm)") #spread is evenish among samples

#the levels of the x variable (soil type) are fixed
#perform Tukey's HSD post-hoc comparisons
TukeyHSD <- HSD.test(lm, "Soil", console=TRUE)
#it's also possible to get the p-values for each comparison
model <- aov(Height~Soil, data=sagedata)
Tukeymodel <- TukeyHSD(model)
Tukeymodel
plot(Tukeymodel) #can visualize the confidence intervals

#prepare data for plotting
graphdata <- sagedata %>%
  group_by(Soil) %>%
  summarize(mean=mean(Height, na.rm=TRUE), se=sd(Height, na.rm=TRUE)/sqrt(length(na.omit(Height))))
graphdata
#rearrange order of the levels on the x-axis
graphdata$Soil <- factor(graphdata$Soil, levels = c("Sterile", "Sand Dune", "Mountains", "Chaparral"))

#make a barplot with ggplot()
ggplot(data=graphdata, aes(x=Soil, y=mean, fill=Soil)) + 
  theme_bw() + 
  theme(axis.text.x=element_text(color="black", size=14), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=16, face="bold"), 
        axis.title.y = element_text(color="black", size=16, face="bold"),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_bar(colour="black", fill="grey", width=0.8, stat="identity") + 
  guides(fill=FALSE) + 
  ylab("Height (cm)") +
  xlab("Soil Type") +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), stat="identity", 
        position=position_dodge(width=0.9), width=0.1) +
  geom_text(aes(label=c("a", "bc", "b", "c"), y = mean + 1))
#note that the labels are assigned in the order of levels in "graphdata"


#Model II ANOVA
#test different coral genets for polyp area

#load data
genets <- read.csv("Data/Genet data.csv")
genets
inspect_df(genets) #check for missing data
#create a univariate ANOVA with genet as a random effect
genets$Genet <- as.factor(genets$Genet)
#lmer() is like lm() but can take random effects
genetmodel <- lmer(Area~1 + (1|Genet), data = genets) #the code 1|Genet sets genet as a random factor
summary(genetmodel) #gives raw variance components

#check assumptions
hist(resid(genetmodel), breaks = 20, xlab = "Area residuals", main = "")
qqPlot(resid(genetmodel))
shapiro.test(resid(genetmodel)) #all okay
#make boxplots of each genet
genet1 <- genets %>%
  filter(Genet == "Genet 1")
genet2 <- genets %>%
  filter(Genet == "Genet 2")
genet3 <- genets %>%
  filter(Genet == "Genet 3")
#groups should have similar dispersion
boxplot(genet1$Area, genet2$Area, genet3$Area, names = c("1", "2", "3"), 
        xlab = "Genet", ylab = "Polyp area") #not exactly perfect
#can also look at a boxplot of pooled sample residuals
boxplot(resid(genetmodel), ylab = "Area residuals") #not too bad
bartlett.test(Area~Genet, data = genets) #not too good
#try a square root transformation
genets$sqrtarea <-sqrt(genets$Area)
genets

#make a glm with sqrt data
sqrtgenetmodel <- lmer(sqrtarea~1 + (1|Genet), data = genets) #the code 1|Genet sets genet as a random factor
summary(sqrtgenetmodel) #gives raw variance components

#check assumptions
hist(resid(sqrtgenetmodel), breaks = 20, xlab = "Sqrt area residuals", main = "")
qqPlot(resid(sqrtgenetmodel))
shapiro.test(resid(sqrtgenetmodel)) #all okay
#make boxplots of each genet
genet1 <- genets %>%
  filter(Genet == "Genet 1")
genet2 <- genets %>%
  filter(Genet == "Genet 2")
genet3 <- genets %>%
  filter(Genet == "Genet 3")
boxplot(genet1$sqrtarea, genet2$sqrtarea, genet3$sqrtarea, names = c("1", "2", "3"), 
        xlab = "Genet", ylab = "Sqrt (polyp area)")
boxplot(resid(sqrtgenetmodel), ylab = "Sqrt area residuals") #can also look at a boxplot of residuals
bartlett.test(sqrtarea~Genet, data = genets) #better but p is still < 0.05

#create a model without genet
genetmodel0 <- lm(sqrtarea ~ 1, data = genets) 
#perform a log likelihood ratio test
anova(sqrtgenetmodel, genetmodel0) #P < 0.05 -- genet has a significant effect
#can also calculate chi-square using the log likelihood scores
X2 <- -2*(89.862-106.995)
1-pchisq(X2, 1) #generate p-value

#calculate variance components
vc <- as.data.frame(VarCorr(sqrtgenetmodel)) 
ICC <- vc$vcov[1] / sum(vc$vcov)
ICC #Intraclass correlation coefficient (%added variance among groups)
varwithin <- 1-ICC #% variance within groups
varwithin

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Area ~ Genet, data = genets) 
#good but really only shows that genet dispersion differs

#perMANOVA
set.seed(123) #makes the permutation output reproducible
adonis2(dist(genets$Area) ~ genets$Genet, permutations = 999) #performs test
#check homogeneity of dispersion
dist_mat <- dist(genets$Area, method = "euclidean") #calculate distance matrix
bd <- betadisper(dist_mat, genets$Genet) #calculate beta dispersion
anova(bd)
permutest(bd)
#tests of homogeneity of dispersion are significant, but perMANOVA isn't detecting an effect of genet
#these perMANOVA results are driven by unequal within-genet variance rather than differences in means
#mixed model ANOVA with sqrt transformed data is the best choice for this analysis


```

---

## assumption_report

```r
choose_column <- function(data, candidates, question) {
  if (length(candidates) == 0) {
    stop("No suitable columns were found for this question.")
  }
  
  choices <- vapply(candidates, function(name) {
    column_number <- match(name, names(data))
    values <- unique(na.omit(data[[name]]))
    examples <- paste(head(format(values, trim = TRUE), 3), collapse = ", ")
    
    paste0(
      "Data column ", column_number, " â€” ", name,
      " [", class(data[[name]])[1], "]",
      if (examples != "") paste0(" | examples: ", examples)
    )
  }, character(1))
  
  selected <- utils::menu(choices, title = question)
  
  if (selected == 0) {
    stop("No column was selected.")
  }
  
  candidates[selected]
}


assumption_report_interactive <- function(data, alpha = 0.05) {
  if (!is.data.frame(data)) {
    stop("Please provide a data frame.")
  }
  
  # Ask which numeric column is the outcome
  numeric_columns <- names(data)[
    vapply(data, is.numeric, logical(1))
  ]
  
  response <- choose_column(
    data,
    numeric_columns,
    "Which numeric column is the outcome?"
  )
  
  # Ask which column identifies the groups
  group_candidates <- names(data)[
    vapply(data, function(x) {
      is.factor(x) || is.character(x) || is.logical(x)
    }, logical(1))
  ]
  
  group <- choose_column(
    data,
    group_candidates,
    "Which column identifies the groups?"
  )
  
  # Ask how observations were collected
  design_choice <- utils::menu(
    c(
      "Independent: different, unrelated subjects in each group",
      "Paired: the same subjects or matched units measured in both groups"
    ),
    title = "How were the observations collected?"
  )
  
  if (design_choice == 0) {
    stop("No study design was selected.")
  }
  
  design <- if (design_choice == 1) "independent" else "paired"
  
  # Keep only the columns needed for this analysis and remove incomplete rows
  d <- data.frame(
    y = data[[response]],
    g = factor(data[[group]])
  )
  
  id <- NULL
  
  if (design == "paired") {
    if (nlevels(droplevels(d$g)) != 2) {
      stop("This paired version requires exactly two groups or conditions.")
    }
    
    id_candidates <- setdiff(names(data), c(response, group))
    
    id <- choose_column(
      data,
      id_candidates,
      "Which column identifies each subject or matched pair?"
    )
    
    d$id <- data[[id]]
  }
  
  d <- d[complete.cases(d), ]
  d$g <- droplevels(d$g)
  
  if (nlevels(d$g) < 2) {
    stop("At least two groups with complete data are needed.")
  }
  
  # Check normality within each independent group,
  # or check the paired differences for paired data
  if (design == "independent") {
    normality_data <- split(d$y, d$g)
  } else {
    group_levels <- levels(d$g)
    ids <- unique(d$id)
    
    if (anyDuplicated(d[c("id", "g")])) {
      stop("There is more than one measurement for a subject in a condition.")
    }
    
    group_1 <- d[d$g == group_levels[1], ]
    group_2 <- d[d$g == group_levels[2], ]
    
    if (!setequal(group_1$id, group_2$id)) {
      stop("The two conditions do not contain the same subjects or pairs.")
    }
    
    values_1 <- group_1$y[match(ids, group_1$id)]
    values_2 <- group_2$y[match(ids, group_2$id)]
    
    normality_data <- list("Paired differences" = values_2 - values_1)
  }
  
  shapiro_p <- vapply(normality_data, function(x) {
    if (length(x) < 3 || length(x) > 5000) {
      return(NA_real_)
    }
    shapiro.test(x)$p.value
  }, numeric(1))
  
  normality_result <- ifelse(
    is.na(shapiro_p),
    "Could not test: sample size must be 3â€“5000",
    ifelse(
      shapiro_p > alpha,
      "No strong evidence against normality",
      "Evidence of departure from normality"
    )
  )
  
  normality_table <- data.frame(
    group_or_measure = names(normality_data),
    n = vapply(normality_data, length, integer(1)),
    shapiro_p = unname(shapiro_p),
    result = unname(normality_result),
    row.names = NULL
  )
  
  # Test roughly equal variances only for independent groups
  variance_p <- NA_real_
  variance_result <- "Not tested for paired data."
  
  if (design == "independent") {
    if (requireNamespace("car", quietly = TRUE)) {
      levene <- car::leveneTest(y ~ g, data = d)
      variance_p <- levene[["Pr(>F)"]][1]
      
      variance_result <- if (variance_p > alpha) {
        "No strong evidence that group variances differ."
      } else {
        "Evidence that group variances differ."
      }
    } else {
      variance_result <- paste(
        "Not tested because the car package is not installed.",
        'To install it, run install.packages("car").'
      )
    }
  }
  
  # Give a next-step suggestion
  if (design == "paired") {
    suggestion <- if (is.na(shapiro_p[1])) {
      "Inspect a Q-Q plot of the paired differences before choosing a test."
    } else if (shapiro_p[1] > alpha) {
      "The checks did not flag a clear problem. Inspect a Q-Q plot of the paired differences, then consider a paired t-test."
    } else {
      "Inspect a Q-Q plot of the paired differences. If the departure is substantial, discuss a Wilcoxon signed-rank test with your instructor."
    }
  } else if (is.na(variance_p)) {
    suggestion <- "Check the group variances and inspect group Q-Q plots before choosing the test."
  } else if (variance_p <= alpha) {
    suggestion <- "Consider Welch's t-test for two groups, or Welch's ANOVA for three or more groups."
  } else {
    suggestion <- "The checks did not flag a clear problem. Inspect group Q-Q plots, then consider the planned independent-samples test."
  }
  
  cat("\nASSUMPTION REPORT\n")
  cat("Outcome column:", response, "\n")
  cat("Group column:", group, "\n")
  cat("Design:", design, "\n\n")
  
  cat("Normality check:\n")
  print(normality_table, row.names = FALSE)
  
  cat("\nEqual-variance check:\n")
  cat(variance_result, "\n")
  if (!is.na(variance_p)) {
    cat("Levene p-value:", signif(variance_p, 4), "\n")
  }
  
  cat("\nIndependence:\n")
  cat("This is based on your answer about how the study was designed.\n")
  
  cat("\nSuggested next step:\n")
  cat(suggestion, "\n")
  
  invisible(list(
    response = response,
    group = group,
    design = design,
    normality = normality_table,
    variance_p = variance_p,
    suggestion = suggestion
  ))
}
```

---

## assumption_report_interactive

```r
# Interactive assumption report for one-sample, independent-group, and paired data.
# Run this file with source(file.choose()) and select this .R file.
# Then run: assumption_report_interactive(mydata)
# Replace mydata with the name of your data frame.
# To choose a CSV without typing its file path, you can use:
# mydata <- read.csv(file.choose())
# If R shows a Browse[1]> prompt after an error, type Q and press Enter to exit it.

assumption_report_interactive <- function(data, alpha = 0.05) {

  # Basic checks before showing any menus
  if (!is.data.frame(data)) {
    stop("The input must be a data frame. Load your CSV first, then pass its name.")
  }
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("The data frame has no rows or columns.")
  }
  if (length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be one number between 0 and 1, such as 0.05.")
  }

  # Print a menu with both the choice number and the original data column number.
  choose_column <- function(candidates, question) {
    if (length(candidates) == 0) return(NULL)

    choices <- vapply(candidates, function(name) {
      column_number <- match(name, names(data))
      values <- unique(na.omit(data[[name]]))
      examples <- paste(head(format(values, trim = TRUE), 3), collapse = ", ")
      type <- paste(class(data[[name]]), collapse = "/")
      paste0(
        "Data column ", column_number, " â€” ", name,
        " [", type, "]",
        if (nzchar(examples)) paste0(" | examples: ", examples)
      )
    }, character(1))

    selected <- utils::menu(choices, title = question)
    if (selected == 0) return(NULL)
    candidates[selected]
  }

  # Flag columns that look like identifiers, so they are not mistaken for outcomes.
  looks_like_id <- function(name, x) {
    name_flag <- grepl("(^id$|individual|subject|sample|tree|plant|specimen)",
                       name, ignore.case = TRUE)
    values <- unique(na.omit(x))
    many_unique <- length(values) >= 0.8 * sum(!is.na(x))
    name_flag || many_unique
  }

  numeric_columns <- names(data)[vapply(data, is.numeric, logical(1))]
  likely_ids <- numeric_columns[vapply(numeric_columns, function(name) {
    looks_like_id(name, data[[name]])
  }, logical(1))]
  outcome_candidates <- setdiff(numeric_columns, likely_ids)

  if (length(outcome_candidates) == 0) {
    outcome_candidates <- numeric_columns
    if (length(outcome_candidates) > 0) {
      cat("All numeric columns look ID-like; review the labels carefully.\n")
    }
  }
  if (length(outcome_candidates) == 0) {
    stop("No numeric columns were found. Check that the outcome was imported as numeric.")
  }

  design_choice <- utils::menu(
    c("One numeric outcome; no groups to compare",
      "Independent groups; different subjects in each group",
      "Paired/repeated data; the same subjects measured in both conditions"),
    title = "What kind of comparison are you doing?"
  )
  if (design_choice == 0) {
    cat("No design selected; no analysis was run.\n")
    return(invisible(NULL))
  }
  design <- c("one-sample", "independent", "paired")[design_choice]

  response <- choose_column(outcome_candidates,
                            "Which numeric column is the outcome?")
  if (is.null(response)) {
    cat("No outcome selected; no analysis was run.\n")
    return(invisible(NULL))
  }

  d <- data.frame(y = data[[response]])
  complete <- !is.na(d$y)
  d <- d[complete, , drop = FALSE]

  group <- NULL
  id <- NULL
  normality_data <- list()
  variance_p <- NA_real_
  variance_text <- "Not applicable for a one-sample or paired comparison."

  if (design == "one-sample") {
    normality_data <- list("Outcome values" = d$y)
  }

  if (design %in% c("independent", "paired")) {
    categorical <- names(data)[vapply(data, function(x) {
      is.factor(x) || is.character(x) || is.logical(x)
    }, logical(1))]
    coded_groups <- names(data)[vapply(data, function(x) {
      is.numeric(x) && length(unique(na.omit(x))) >= 2 &&
        length(unique(na.omit(x))) <= 10
    }, logical(1))]
    group_candidates <- setdiff(unique(c(categorical, coded_groups)), response)
    group_candidates <- group_candidates[!vapply(group_candidates, function(name) {
      looks_like_id(name, data[[name]])
    }, logical(1))]

    if (length(group_candidates) == 0) {
      cat("\nI couldn't find a suitable group column in this data.\n")
      cat("This often means the data are one sample, or the group column needs to be converted to a factor.\n")
      cat("For example, if a numeric column contains group codes, run: data$Group <- factor(data$Group)\n")
      cat("No report was generated for a group comparison.\n")
      return(invisible(list(
        response = response,
        design = design,
        message = "No suitable group column found."
      )))
    }

    group <- choose_column(group_candidates,
                           "Which column identifies the groups or conditions?")
    if (is.null(group)) {
      cat("No group selected; no analysis was run.\n")
      return(invisible(NULL))
    }

    d$g <- factor(data[[group]][complete])
    d <- d[!is.na(d$g), , drop = FALSE]
    d$g <- droplevels(d$g)

    if (nlevels(d$g) < 2) {
      cat("The selected group column has fewer than two groups after removing missing rows.\n")
      return(invisible(NULL))
    }

    if (design == "independent") {
      normality_data <- split(d$y, d$g)
      if (all(table(d$g) >= 2)) {
        if (requireNamespace("car", quietly = TRUE)) {
          levene <- tryCatch(car::leveneTest(y ~ g, data = d), error = function(e) NULL)
          if (!is.null(levene)) {
            variance_p <- levene[["Pr(>F)"]][1]
            variance_text <- if (is.finite(variance_p) && variance_p > alpha) {
              "No strong evidence that group variances differ."
            } else if (is.finite(variance_p)) {
              "Evidence that group variances differ."
            } else {
              "Levene's test did not return a usable p-value. Check group sizes and data."
            }
          } else {
            variance_text <- "Levene's test could not run; check group sizes and data."
          }
        } else {
          variance_text <- paste(
            "Not tested: install the car package to run Levene's test.",
            'Use install.packages("car") once, then rerun the function.'
          )
        }
      } else {
        variance_text <- "Not tested: each group needs at least two observations for this check."
      }
    }

    if (design == "paired") {
      if (nlevels(d$g) != 2) {
        cat("Paired comparisons require exactly two conditions in this version.\n")
        return(invisible(NULL))
      }

      id_candidates <- setdiff(names(data), c(response, group))
      if (length(id_candidates) == 0) {
        cat("No possible subject/pair ID column was found.\n")
        return(invisible(NULL))
      }
      id <- choose_column(id_candidates,
                          "Which column identifies the same subject or pair?")
      if (is.null(id)) {
        cat("No ID column selected; no analysis was run.\n")
        return(invisible(NULL))
      }

      d$id <- data[[id]][complete][!is.na(data[[group]][complete])]
      d <- d[!is.na(d$id), , drop = FALSE]
      if (nrow(d) == 0) {
        cat("No complete subject, condition, and outcome rows remain.\n")
        return(invisible(NULL))
      }
      if (anyDuplicated(d[c("id", "g")])) {
        cat("There is more than one outcome per subject per condition. Check the data layout.\n")
        return(invisible(NULL))
      }

      group_levels <- levels(d$g)
      first <- d[d$g == group_levels[1], , drop = FALSE]
      second <- d[d$g == group_levels[2], , drop = FALSE]
      if (!setequal(first$id, second$id)) {
        cat("The two conditions do not contain the same subjects/pairs.\n")
        return(invisible(NULL))
      }
      paired_differences <- second$y[match(first$id, second$id)] - first$y
      normality_data <- list("Paired differences" = paired_differences)
    }
  }

  # Shapiro-Wilk can only be run for samples with 3 to 5000 observations.
  sample_sizes <- vapply(normality_data, length, integer(1))
  shapiro_p <- vapply(normality_data, function(x) {
    if (length(x) < 3 || length(x) > 5000) return(NA_real_)
    tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
  }, numeric(1))

  normality_result <- ifelse(
    is.na(shapiro_p),
    "Could not test (Shapiro-Wilk needs n = 3 to 5000)",
    ifelse(shapiro_p > alpha,
           "No strong evidence against normality",
           "Evidence of departure from normality")
  )

  normality_table <- data.frame(
    data_checked = names(normality_data),
    n = sample_sizes,
    shapiro_p = unname(shapiro_p),
    result = unname(normality_result),
    row.names = NULL
  )

  if (design == "one-sample") {
    suggestion <- paste(
      "Check a histogram or Q-Q plot and confirm observations are independent.",
      "If the outcome is reasonably normal and a one-sample mean matches your question,",
      "a one-sample t-test may be appropriate."
    )
  } else if (design == "paired") {
    suggestion <- if (!is.na(shapiro_p[1]) && shapiro_p[1] <= alpha) {
      "Inspect a Q-Q plot of the paired differences. If departure is substantial, ask whether a paired nonparametric method is appropriate."
    } else {
      "Inspect a Q-Q plot of the paired differences. If reasonable, a paired t-test may be appropriate."
    }
  } else if (is.na(variance_p)) {
    suggestion <- "Inspect group Q-Q plots and variances. Consider Welch's test if equal variances are doubtful."
  } else if (variance_p <= alpha) {
    suggestion <- if (nlevels(d$g) == 2) {
      "Consider Welch's independent-samples t-test, which does not assume equal variances."
    } else {
      "Consider Welch's ANOVA, which does not assume equal variances."
    }
  } else {
    suggestion <- "Inspect group Q-Q plots and confirm independence from the study design before using the planned independent-groups test."
  }

  cat("\nASSUMPTION REPORT\n")
  cat("Outcome:", response, "(data column", match(response, names(data)), ")\n")
  if (!is.null(group)) {
    cat("Groups/conditions:", group, "(data column", match(group, names(data)), ")\n")
  }
  if (!is.null(id)) {
    cat("Pair ID:", id, "(data column", match(id, names(data)), ")\n")
  }
  cat("Design:", design, "\n\n")
  cat("Normality check:\n")
  print(normality_table, row.names = FALSE)
  cat("\nEqual-variance check:\n", variance_text, "\n", sep = "")
  if (!is.na(variance_p)) cat("Levene p-value:", signif(variance_p, 4), "\n")
  cat("\nIndependence:\n")
  cat("This cannot be confirmed from the values alone; verify it from how the data were collected.\n")
  cat("\nSuggested next step:\n", suggestion, "\n", sep = "")
  cat("\nA p-value above alpha means the test did not find strong evidence against an assumption; it does not prove the assumption is true.\n")

  invisible(list(
    response = response,
    group = group,
    id = id,
    design = design,
    normality = normality_table,
    variance_p = variance_p,
    suggestion = suggestion
  ))
}

```

---

## Assumptions Testing

```r
assumption_report <- function(data, response, group,
                              alpha = 0.05,
                              independence = "unknown") {
  #assumption_report(data = mydata, response = "Height",group = "Strain", independence = "yes" )
  # Pull out the outcome and group columns, and remove incomplete rows
  d <- data.frame(
    y = data[[response]],
    g = factor(data[[group]])
  )
  
  d <- d[complete.cases(d), ]
  d$g <- droplevels(d$g)
  
  if (!is.numeric(d$y)) {
    stop("The response column must be numeric.")
  }
  
  if (nlevels(d$g) < 2) {
    stop("The group column must contain at least two groups.")
  }
  
  if (!independence %in% c("yes", "no", "unknown")) {
    stop('Use independence = "yes", "no", or "unknown".')
  }
  
  # Shapiro-Wilk normality test within each group
  group_data <- split(d$y, d$g)
  
  normality_p <- sapply(group_data, function(x) {
    if (length(x) < 3 || length(x) > 5000) {
      return(NA_real_)
    }
    
    shapiro.test(x)$p.value
  })
  
  normality_result <- ifelse(
    is.na(normality_p),
    "Could not test: group size must be 3â€“5000",
    ifelse(
      normality_p > alpha,
      "No strong evidence against normality",
      "Evidence of departure from normality"
    )
  )
  
  normality_table <- data.frame(
    group = names(group_data),
    n = sapply(group_data, length),
    shapiro_p = as.numeric(normality_p),
    result = normality_result,
    row.names = NULL
  )
  
  # Levene's test for roughly equal group variances
  if (requireNamespace("car", quietly = TRUE)) {
    levene_result <- car::leveneTest(y ~ g, data = d)
    variance_p <- levene_result[["Pr(>F)"]][1]
    
    variance_text <- if (variance_p > alpha) {
      "No strong evidence that the group variances differ."
    } else {
      "Evidence that the group variances differ."
    }
  } else {
    levene_result <- NULL
    variance_p <- NA_real_
    variance_text <- paste(
      "Not tested. Install the car package to run Levene's test:",
      'install.packages("car")'
    )
  }
  
  # Suggested next step
  if (independence == "no") {
    suggestion <- paste(
      "The observations are not independent. Recheck the study design;",
      "a paired or repeated-measures analysis may be needed."
    )
  } else if (variance_p <= alpha) {
    suggestion <- paste(
      "Consider a method that does not assume equal variances,",
      "such as Welch's t-test for two independent groups."
    )
  } else if (independence == "unknown") {
    suggestion <- paste(
      "Confirm independence from the study design.",
      "Also inspect histograms or Q-Q plots for strong skew or outliers."
    )
  } else {
    suggestion <- paste(
      "These checks did not flag a clear problem.",
      "Inspect histograms or Q-Q plots, then proceed with the planned test"
    )
  }
  
  report <- c(
    "ASSUMPTION REPORT",
    paste("Response:", response),
    paste("Groups:", group),
    paste("Rows analyzed:", nrow(d)),
    "",
    "Normality by group (Shapiro-Wilk):",
    capture.output(print(normality_table, row.names = FALSE)),
    "",
    paste("Equal variances (Levene's test, p =",
          ifelse(is.na(variance_p), "not available",
                 signif(variance_p, 3)), "):"),
    variance_text,
    "",
    paste("Independence:", independence),
    if (independence == "unknown") {
      "Check this using how the observations were collected."
    } else {
      ""
    },
    "",
    "Suggested next step:",
    suggestion
  )
  
  cat(paste(report, collapse = "\n"))
  
  invisible(list(
    report = report,
    normality = normality_table,
    levene = levene_result,
    variance_p = variance_p
  ))
}
```

---

## Biostats_330_Master_Analysis_Workflow

```r
# Biostats 330 â€” Master Analysis Workflow
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

# NOTE â€” What is the first thing to do?
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

# NOTE â€” How should I import the data?
# Easiest answer: use file.choose() when you are unsure about the path.
# For a reproducible submission, replace it with a project-relative path such
# as "data/Orange trees FA26 - Sheet1.csv".

data <- read.csv(file.choose(), check.names = FALSE)

# For Excel files, use this instead:
# data <- readxl::read_excel(file.choose())

# NOTE â€” What should I look at immediately after importing?
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

# NOTE â€” What if the file will not load?
#   - If the error says the file cannot be found, use file.choose().
#   - If the columns are in one giant column, check the delimiter and use
#     read.csv2() or read.delim() only if the file actually uses that format.
#   - If a column name contains spaces, use data[["Column name"]] or rename it.
#   - In this repository, use lowercase data/ rather than Data/ when possible.


# ============================================================================
# 2. WRITE THE QUESTION BEFORE RUNNING A TEST
# ============================================================================

# NOTE â€” What am I trying to compare or relate?
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

# NOTE â€” How do I choose "less", "greater", or "two.sided"?
# Ask: Did the question predict a direction before seeing the data?
#   - Yes, lower -> alternative = "less"
#   - Yes, higher -> alternative = "greater"
#   - No specific direction, or the professor says two-tailed -> "two.sided"
# Do not change the direction after seeing the result.


# ============================================================================
# 3. IDENTIFY THE DESIGN â€” THE DECISION THAT CONTROLS EVERYTHING AFTER IT
# ============================================================================

# NOTE â€” Are the observations independent?
# Ask: Could one observation's measurement influence another observation's
# measurement, or were both measurements taken from the same subject/tree/deer?
#   - Same individual or deliberately matched units -> paired/repeated design.
#   - Different subjects/trees/plots with no matching -> independent design.
#   - If unsure, look at the study description, not at the p-value.
# Independence is determined by how the data were collected; R cannot repair
# a design that is not independent.

# NOTE â€” Is the response numeric or categorical?
#   - Numeric: length, mass, height, temperature, concentration, diversity.
#   - Categorical/count: red/black, rural/industrial, observed category counts.
#   - ID columns such as Individual, Subject, Tree number, or Quadrat are not
#     outcomes. Do not test an ID column.


# ============================================================================
# 4. DESCRIPTIVE GRAPHS BEFORE FORMAL ASSUMPTION TESTS
# ============================================================================

# NOTE â€” What should I graph?
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
# 5. ASSUMPTIONS â€” SIMPLE QUESTIONS AND WHAT EACH ANSWER MEANS
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

# NOTE â€” What does Shapiro-Wilk test?
# H_0: the values are consistent with a normal distribution.
# If p <= 0.05, there is evidence against normality.
# If p > 0.05, you did not detect evidence against normality; this does not
# prove the data are perfectly normal. Always inspect the histogram and Q-Q
# plot too, especially with small samples.

# NOTE â€” Which values should I test for normality?
#   - One-sample test: the response values.
#   - Independent two-group test: the response within each group.
#   - Paired test: the paired differences, not the two raw columns separately.
#   - ANOVA: the model residuals, not only the raw response.
#   - Regression: the model residuals, plus the residual-vs-fitted plot.
#   - Repeated measures ANOVA: the model residuals and the study design.

# Example:
# normality_report(data$RESPONSE, "RESPONSE")
# normality_report(data$RESPONSE[data$GROUP == "Group 1"], "Group 1")

# NOTE â€” How do I check equal variances?
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

# NOTE â€” Do I need a transformation?
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

# NOTE â€” Which version should I keep?
# Keep the version with the most reasonable biological scale and the best
# overall diagnostics: Q-Q plot, histogram, residual pattern, spread, and
# skewness. Shapiro-Wilk alone does not choose the transformation.
# The final statistical test must use the chosen version consistently.


# ============================================================================
# 7. BRANCH A â€” ONE NUMERIC SAMPLE AGAINST A FIXED VALUE
# ============================================================================

# Used in your work for questions such as Buckwheat density compared with 4.

# NOTE â€” Is the question "is the population mean different from VALUE"?
#   - Yes -> one-sample t-test if the response is reasonably normal.
#   - If the response remains strongly nonnormal, ask whether the professor
#     wants a one-sample sign/Wilcoxon method.

# H_0: mean(RESPONSE) = VALUE
# H_A: mean(RESPONSE) != VALUE  # change direction only if predicted

# t.test(data$RESPONSE, mu = VALUE, alternative = "two.sided")

# Report: sample mean, confidence interval, t, df, p-value, and the decision
# about H_0. If significant, state the direction using the sample mean.


# ============================================================================
# 8. BRANCH B â€” ONE PROPORTION OR ONE CATEGORICAL RATIO
# ============================================================================

# B1. Exact binomial test â€” used for the Orange tree fruit question.

# NOTE â€” Is each observation a success/failure and is the null a probability?
#   - Yes -> use binom.test().
#   - Count successes with sum(), count trials with nrow() or the stated n.

# total <- nrow(data)
# successes <- sum(data$SUCCESS_COLUMN == 1, na.rm = TRUE)
# binom.test(successes, total, p = 0.95, alternative = "less")

# NOTE â€” What does the alternative mean here?
# "less" asks whether the true success rate is below the stated rate.
# "greater" asks whether it is above the stated rate.
# "two.sided" asks whether it differs in either direction.

# B2. Goodness-of-fit chi-square/G-test â€” used for bee sex ratio and the
# dihybrid 9:3:3:1 question.

# NOTE â€” Is one categorical variable being compared with a theoretical ratio?
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

# NOTE â€” What do I conclude?
#   - p <= alpha: reject H_0; the observed counts do not fit the proposed
#     ratio closely enough.
#   - p > alpha: fail to reject H_0; there is not enough evidence that the
#     observed counts differ from the proposed ratio.
# Do not say "the null is proven true."

# B3. Intrinsic/Poisson goodness-of-fit â€” used for count frequencies such as
# oranges per tree and the Midterm Q3/Mites-and-midges style problem.

# NOTE â€” Is the question whether count frequencies follow a Poisson model?
#   - Yes -> estimate lambda from the observed counts, calculate expected
#     Poisson probabilities, convert them to expected counts, and combine rare
#     high-count categories (for example, a final "4+" category).
#   - When lambda is estimated from the same data, use df = final number of
#     categories - 2, as in the professor's notes.
#   - Then compare the observed and expected counts with the required G-test
#     or chi-square calculation.


# ============================================================================
# 9. BRANCH C â€” TWO CATEGORICAL VARIABLES: CONTINGENCY TABLES
# ============================================================================

# Used for Ladybird color by location, Foraging behavior by substrate, and
# bird-guild style data.

# NOTE â€” Does the question ask whether two categories are associated?
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

# NOTE â€” Are the expected counts large enough for the planned chi-square/G-test?
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

# NOTE â€” What should the graph show?
#   - Counts: side-by-side bars for raw observed counts.
#   - Proportions: stacked/fill bars when comparing composition within rows.
# The graph describes the pattern; the test decides whether the association
# is stronger than expected by chance.


# ============================================================================
# 10. BRANCH D â€” TWO INDEPENDENT NUMERIC GROUPS
# ============================================================================

# Used in your work for examples such as Microalgae concentration and the
# Barnacle-style two-sample script.

# NOTE â€” Are the two groups independent?
#   - Yes -> continue here.
#   - No, same individuals/paired units -> go to Branch E.

# First inspect each group:
# boxplot(RESPONSE ~ GROUP, data = data)
# normality_report(data$RESPONSE[data$GROUP == "GROUP_1"], "GROUP_1")
# normality_report(data$RESPONSE[data$GROUP == "GROUP_2"], "GROUP_2")
# car::leveneTest(RESPONSE ~ GROUP, data = data)

# NOTE â€” Which test do I run?
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

# NOTE â€” Why does var.equal matter?
# var.equal = TRUE assumes the two population variances are equal.
# var.equal = FALSE is Welch's method and is safer when the variances differ.
# It does not mean the data are nonparametric; it is still a t-test for means.


# ============================================================================
# 11. BRANCH E â€” TWO PAIRED NUMERIC MEASUREMENTS
# ============================================================================

# Used for deer hindleg vs foreleg, butterfly wings, and paired cone/condition
# examples. The pairing must be preserved by row.

# NOTE â€” Is row 1's first measurement matched to row 1's second measurement?
#   - Yes -> paired analysis.
#   - No -> independent-groups analysis instead.

# differences <- data$MEASUREMENT_1 - data$MEASUREMENT_2
## normality_report(differences, "paired differences")

# NOTE â€” What is the assumption decision based on?
# The normality of the DIFFERENCES, not whether each raw column is normal.

# If differences are reasonably normal:
# t.test(data$MEASUREMENT_1, data$MEASUREMENT_2,
#        paired = TRUE, alternative = "two.sided")

# If differences remain strongly nonnormal:
# wilcox.test(data$MEASUREMENT_1, data$MEASUREMENT_2,
#             paired = TRUE, alternative = "two.sided")

# If the professor specifically asks for a paired sign test, use the sign of
# each nonzero difference and report that method explicitly.

# NOTE â€” What if the assignment gives an explicit direction?
# Follow the assignment's direction. For example, PSet2Q5 explicitly called
# for a two-tailed paired t-test even though the biological wording suggested
# a directional expectation.


# ============================================================================
# 12. BRANCH F â€” THREE OR MORE INDEPENDENT NUMERIC GROUPS
# ============================================================================

# Used for Erodium, Aphid length by gall, diversity by location, Sagebrush,
# and similar ANOVA questions.

# NOTE â€” Are there three or more independent groups?
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

# NOTE â€” What is the simple decision?
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

# NOTE â€” What if there is a random factor?
# Your ANOVA script includes Genet as a random factor. Do not replace that
# design with an ordinary one-way ANOVA without checking the study structure.
# The relevant model pattern is:
# lme4::lmer(RESPONSE ~ 1 + (1 | RANDOM_FACTOR), data = data)


# ============================================================================
# 13. BRANCH G â€” REPEATED-MEASURES ANOVA
# ============================================================================

# Used for Stress/BPM and Coral/FvFm questions.

# NOTE â€” Did the same subject/coral/tree receive multiple conditions or days?
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

# NOTE â€” What does the ID column do?
# It tells R which rows belong to the same individual. Without it, R cannot
# model the repeated measurements correctly.

# NOTE â€” What comes after a significant repeated-measures ANOVA?
# Use emmeans to compare the conditions and explain which condition means
# differ. Plot the means with uncertainty bars, but do not use the graph alone
# to decide significance.


# ============================================================================
# 14. BRANCH H â€” CORRELATION AND LINEAR REGRESSION
# ============================================================================

# Used for penguin beak width/depth, environmental variables, ants, transects,
# and stream biomass/species richness.

# NOTE â€” Is the question about association between two numeric variables?
#   - Yes -> make a scatterplot first.

# plot(data$X, data$Y, pch = 19, xlab = "X", ylab = "Y")

# NOTE â€” Does the scatterplot look approximately linear, with no extreme
# outlier controlling the pattern?
#   - Yes -> Pearson correlation is reasonable.
#   - No, but the pattern is monotonic -> Spearman correlation may be better.
#   - No clear pattern -> do not force a correlation just because both columns
#     are numeric.

# Pearson:
# cor.test(data$X, data$Y, method = "pearson")

# Spearman:
# cor.test(data$X, data$Y, method = "spearman")

# NOTE â€” Do I need a linear model too?
#   - If you only need strength/direction of association -> cor.test may be
#     enough.
#   - If you need a fitted equation, slope, prediction, or residual diagnostics
#     -> fit lm().

# model <- lm(Y ~ X, data = data)
# summary(model)
# plot(model, which = 1)  # residuals vs fitted: random cloud around 0
# plot(model, which = 2)  # normal Q-Q plot of residuals
# shapiro.test(resid(model))

# NOTE â€” What do I look for in the residual-vs-fitted plot?
#   - Random cloud around zero -> linearity/spread look reasonable.
#   - Curved pattern -> relationship may not be linear.
#   - Funnel shape -> spread changes with fitted values.
#   - One point far away -> possible influential outlier.

# For many correlations, as in PSet6Q2, adjust the p-values:
# variables <- c("X1", "X2", "X3", "X4")
# pairs <- combn(variables, 2, simplify = FALSE)
# results$p_adjusted <- p.adjust(results$p, method = "holm")


# ============================================================================
# 15. WIDE VERSUS LONG DATA â€” ONLY RESHAPE WHEN THE METHOD NEEDS IT
# ============================================================================

# NOTE â€” Is one row already one subject and separate columns are conditions?
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

# NOTE â€” What should stay the same after reshaping?
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
# #notes â€” p-value decision:
# "Because p = [p] is [less/greater] than alpha = 0.05, I [reject/fail to
# reject] H_0."
#
# #notes â€” numeric comparison:
# "The [test] showed that [response] was [higher/lower/different] in [group]
# than [comparison], [statistic] = [value], df = [value], p = [value]."
#
# #notes â€” correlation:
# "There was a [positive/negative] [weak/moderate/strong] association between
# X and Y, r = [r], p = [p]. This describes association, not causation."
#
# #notes â€” non-significant result:
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


```

---

## Contingency tables practice extra snippets

```r


# View observed contingency table
contable_obs

# Add across each ROW
rowSums(contable_obs)

# Add down each COLUMN
colSums(contable_obs)

# Add EVERY observation = grand total
sum(contable_obs)

# Display contingency table WITH row and column totals
addmargins(contable_obs)

#calculate totals
redtotal <- sum(redrural, redindustrial)
blacktotal <- sum(blackrural, blackindustrial)
ruraltotal <- sum(redrural, blackrural)
industrialtotal <- sum(redindustrial, blackindustrial)
rowcoltotal <- redtotal+blacktotal

#make a contingency table of observed values
contable_obs <- matrix(c(redrural, redindustrial, redtotal, blackrural, blackindustrial, blacktotal,
                         ruraltotal, industrialtotal, rowcoltotal), nrow=3, ncol=3)
contable_obs
#name rows and columns
row.names(contable_obs) <- (c("Rural", "Industrial", "Total"))
colnames(contable_obs) <- (c("Red", "Black", "Total"))
contable_obs

# Add across each ROW
rowSums(contable_obs)

# Add down each COLUMN
colSums(contable_obs)

# Add every observation in the table
sum(contable_obs)

addmargins(contable_obs)

#calculate expected values
redrural_exp <- (100/300)*(155/300)*300
redindustrial_exp <- (145/300)*(155/300)*300
blackrural_exp <- (100/300)*(145/300)*300
blackindustrial_exp <- (145/300)*(145/300)*300

#make a contingency table of expected values
contable_exp <- matrix(c(redrural_exp, redindustrial_exp, redtotal, blackrural_exp, blackindustrial_exp, 
                         blacktotal, ruraltotal, industrialtotal, rowcoltotal), nrow=3, ncol=3)
contable_exp
#name rows and columns
row.names(contable_exp) <- (c("Rural", "Industrial", "Total"))
colnames(contable_exp) <- (c("Red", "Black", "Total"))
contable_exp

#test of independence
#reassemble contingency table of observed values without totals
contable <- matrix(c(redrural, redindustrial, blackrural, blackindustrial), nrow=2, ncol=2)
row.names(contable) <- (c("Rural", "Industrial"))
colnames(contable) <- (c("Red", "Black"))
#perform G test
GTest(contable)

#make plot
names <- matrix(c("Rural", "Ind"))
barplot(contable, beside = TRUE, names.arg = names(contable), legend.text = (names), 
        bty = "n", ylim = c(0, 40), ylab = "Counts", args.legend = list(x = "topright", 
                                                                        inset = c(-0.1, -0.2),bty = "n", cex = 1))
```

---

## contingency tables2

```r
#Contingency tables
#Imports Spreadsheet
data <- read.csv("Data/Ladybird data.csv")

#Opens Spreadsheet In New Tab
View(data)

#load libraries 
library(tidyverse)
library(DescTools) #for G test
library(ggplot2) #for plotting

#analyze dispersal and reproduction data

#create data
redrural <- 70
redindustrial <- 85
blackrural <- 30
blackindustrial <- 115

#calculate totals
redtotal <- sum(redrural, redindustrial)
blacktotal <- sum(blackrural, blackindustrial)
ruraltotal <- sum(redrural, blackrural)
industrialtotal <- sum(redindustrial, blackindustrial)
rowcoltotal <- redtotal+blacktotal

#make a contingency table of observed values
contable_obs <- matrix(c(redrural, redindustrial, redtotal, blackrural, blackindustrial, blacktotal,
                   ruraltotal, industrialtotal, rowcoltotal), nrow=3, ncol=3)
contable_obs
#name rows and columns
row.names(contable_obs) <- (c("Rural", "Industrial", "Total"))
colnames(contable_obs) <- (c("Red", "Black", "Total"))
contable_obs

#calculate expected values
redrural_exp <- (100/300)*(155/300)*300
redindustrial_exp <- (145/300)*(155/300)*300
blackrural_exp <- (100/300)*(145/300)*300
blackindustrial_exp <- (145/300)*(145/300)*300

#make a contingency table of expected values
contable_exp <- matrix(c(redrural_exp, redindustrial_exp, redtotal, blackrural_exp, blackindustrial_exp, 
                       blacktotal, ruraltotal, industrialtotal, rowcoltotal), nrow=3, ncol=3)
contable_exp
#name rows and columns
row.names(contable_exp) <- (c("Rural", "Industrial", "Total"))
colnames(contable_exp) <- (c("Red", "Black", "Total"))
contable_exp

#test of independence
#reassemble contingency table of observed values without totals
contable <- matrix(c(redrural, redindustrial, blackrural, blackindustrial), nrow=2, ncol=2)
row.names(contable) <- (c("Rural", "Industrial"))
colnames(contable) <- (c("Red", "Black"))
#perform G test
GTest(contable)

#make plot
names <- matrix(c("Rural", "Ind"))
barplot(contable, beside = TRUE, names.arg = names(contable), legend.text = (names), 
        bty = "n", ylim = c(0, 40), ylab = "Counts", args.legend = list(x = "topright", 
        inset = c(-0.1, -0.2),bty = "n", cex = 1))
```

---

## Correlation

```r
#Correlation

#clear the environment
rm(list=ls())

#load libraries
library(tidyverse)
library(car)
library(ggplot2)
library(GGally)
library(psych)

#analyze desert ant nest distributions 
#test for correlation between mean distance to vegetation and neighbor distances

#load data
antdata <- read.csv("Data/Ant data.csv")
View(antdata)
str(antdata) #look at structure of variables in data frame

#there's alot of info here
#let's work with the columns "Mean" and "Neighbors"

#check assumptions of normality
#start with mean data
hist(antdata$Mean)
shapiro.test(antdata$Mean)
qqPlot(antdata$Mean)
#the data for mean distance is not normally distributed

#try a log transformation
antdata$logmean <- log(antdata$Mean)
hist(antdata$logmean)
shapiro.test(antdata$logmean)
qqPlot(antdata$logmean)
#log transformation improves normality

#what about square root transformation
antdata$sqrtmean <- sqrt(antdata$Mean)
hist(antdata$sqrtmean)
shapiro.test(antdata$sqrtmean)
qqPlot(antdata$sqrtmean)
#not as good as the log transformation

#how about neighbor distance data
hist(antdata$Neighbors)
shapiro.test(antdata$Neighbors) 
qqPlot(antdata$Neighbors)
#also pretty bad

#let's try a log transformation
antdata$logneighbors <- log(antdata$Neighbors)
hist(antdata$logneighbors)
shapiro.test(antdata$logneighbors) 
qqPlot(antdata$logneighbors)
#log transformation improves normality

#what about square root transformation
antdata$sqrtneighbors <- sqrt(antdata$Neighbors)
hist(antdata$sqrtneighbors)
shapiro.test(antdata$sqrtneighbors)
qqPlot(antdata$sqrtneighbors)
#not as good as the log transformation

#let's try Pearson's correlation with and without the log transformed data
cor.test(antdata$logneighbors, antdata$logmean, method="pearson", na.rm=TRUE)
cor.test(antdata$Neighbors, antdata$Mean, method="pearson", na.rm=TRUE)
#it's a highly significant association either way
#the log transformation has slightly better model fit (r = 0.80 vs. 0.77)

#how about Spearman's rank test (non-parametric) with raw data
cor.test(antdata$Neighbors, antdata$Mean, method="spearman", na.rm=TRUE)
#it gives a warning about ties (the ranks don't fit perfectly)
#interpret the P value as an approximation

#make a scatterplot
ggplot(antdata, aes(x=logmean, y=logneighbors)) +
  geom_point(size = 2) +
  xlab("Log (Mean distance to vegetation)") +
  ylab("Log (Mean distance to neighbors)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
#we could color code the points by site
ggplot(antdata, aes(x=logmean, y=logneighbors, col = Site)) +
  geom_point(size = 2) +
  xlab("Log (Mean distance to vegetation)") +
  ylab("Log (Mean distance to neighbors)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())

#let's make a SPLOM with Transverse range transect data
#we'll look at shrub diversity, evenness, percent cover, and species richness

#load data
transects <- read.csv("Data/Transverse transect data.csv")
View(transects)
str(transects) #look at structure of variables in data frame

#we should probably check the model assumptions
#let's skip that here as the code to do so is available above

#select columns
data_subset <- transects[, c("Diversity", "Evenness", "Richness", "Cover")]

#create the scatter plot matrix
pairs(data_subset, 
      main = "SPLOM",
      pch = 20,                                     
      cex = 1.5,                                    
      col = "#69b3a2"                               
)

#alternate method using ggpairs() from library(GGgally) that looks a little better
#the diagonal plots are histograms
ggpairs(data_subset, 
        title = "")

#yet another method using pair.panels() from library(psych)
#this gives histograms and ellipses
pairs.panels(data_subset,
             method = "pearson", # Correlation method (e.g., pearson, spearman)
             hist.col = "#00AFBB", # Color for histograms
             density = TRUE,      # Add density plots
             ellipses = TRUE      # Add correlation ellipses
)

#create the individual Pearson's correlation tests
#notice that the results match those shown above
cor.test(transects$Diversity, transects$Evenness, method="pearson", na.rm=TRUE)
cor.test(transects$Diversity, transects$Richness, method="pearson", na.rm=TRUE)
cor.test(transects$Diversity, transects$Cover, method="pearson", na.rm=TRUE)
cor.test(transects$Evenness, transects$Richness, method="pearson", na.rm=TRUE)
cor.test(transects$Evenness, transects$Cover, method="pearson", na.rm=TRUE)
cor.test(transects$Cover, transects$Richness, method="pearson", na.rm=TRUE)

```

---

## correlationn

```r
#Correlation

#clear the environment
rm(list=ls())

#load libraries
library(tidyverse)
library(car)
library(ggplot2)
library(GGally)
library(psych)

#analyze desert ant nest distributions 
#test for correlation between mean distance to vegetation and neighbor distances

#load data
antdata <- read.csv("Data/Ant data.csv")
View(antdata)
str(antdata) #look at structure of variables in data frame

#there's alot of info here
#let's work with the columns "Mean" and "Neighbors"

#check assumptions of normality
#start with mean data
hist(antdata$Mean)
shapiro.test(antdata$Mean)
qqPlot(antdata$Mean)
#the data for mean distance is not normally distributed

#try a log transformation
antdata$logmean <- log(antdata$Mean)
hist(antdata$logmean)
shapiro.test(antdata$logmean)
qqPlot(antdata$logmean)
#log transformation improves normality

#what about square root transformation
antdata$sqrtmean <- sqrt(antdata$Mean)
hist(antdata$sqrtmean)
shapiro.test(antdata$sqrtmean)
qqPlot(antdata$sqrtmean)
#not as good as the log transformation

#how about neighbor distance data
hist(antdata$Neighbors)
shapiro.test(antdata$Neighbors) 
qqPlot(antdata$Neighbors)
#also pretty bad

#let's try a log transformation
antdata$logneighbors <- log(antdata$Neighbors)
hist(antdata$logneighbors)
shapiro.test(antdata$logneighbors) 
qqPlot(antdata$logneighbors)
#log transformation improves normality

#what about square root transformation
antdata$sqrtneighbors <- sqrt(antdata$Neighbors)
hist(antdata$sqrtneighbors)
shapiro.test(antdata$sqrtneighbors)
qqPlot(antdata$sqrtneighbors)
#not as good as the log transformation

#let's try Pearson's correlation with and without the log transformed data
cor.test(antdata$logneighbors, antdata$logmean, method="pearson", na.rm=TRUE)
cor.test(antdata$Neighbors, antdata$Mean, method="pearson", na.rm=TRUE)
#it's a highly significant association either way
#the log transformation has slightly better model fit (r = 0.80 vs. 0.77)

#how about Spearman's rank test (non-parametric) with raw data
cor.test(antdata$Neighbors, antdata$Mean, method="spearman", na.rm=TRUE)
#it gives a warning about ties (the ranks don't fit perfectly)
#interpret the P value as an approximation

#make a scatterplot
ggplot(antdata, aes(x=logmean, y=logneighbors)) +
  geom_point(size = 2) +
  xlab("Log (Mean distance to vegetation)") +
  ylab("Log (Mean distance to neighbors)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
#we could color code the points by site
ggplot(antdata, aes(x=logmean, y=logneighbors, col = Site)) +
  geom_point(size = 2) +
  xlab("Log (Mean distance to vegetation)") +
  ylab("Log (Mean distance to neighbors)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())

#let's make a SPLOM with Transverse range transect data
#we'll look at shrub diversity, evenness, percent cover, and species richness

#load data
transects <- read.csv("Data/Transverse transect data.csv")
View(transects)
str(transects) #look at structure of variables in data frame

#we should probably check the model assumptions
#let's skip that here as the code to do so is available above

#select columns
data_subset <- transects[, c("Diversity", "Evenness", "Richness", "Cover")]

#create the scatter plot matrix
pairs(data_subset, 
      main = "SPLOM",
      pch = 20,                                     
      cex = 1.5,                                    
      col = "#69b3a2"                               
)

#alternate method using ggpairs() from library(GGgally) that looks a little better
#the diagonal plots are histograms
ggpairs(data_subset, 
        title = "")

#yet another method using pair.panels() from library(psych)
#this gives histograms and ellipses
pairs.panels(data_subset,
             method = "pearson", # Correlation method (e.g., pearson, spearman)
             hist.col = "#00AFBB", # Color for histograms
             density = TRUE,      # Add density plots
             ellipses = TRUE      # Add correlation ellipses
)

#create the individual Pearson's correlation tests
#notice that the results match those shown above
cor.test(transects$Diversity, transects$Evenness, method="pearson", na.rm=TRUE)
cor.test(transects$Diversity, transects$Richness, method="pearson", na.rm=TRUE)
cor.test(transects$Diversity, transects$Cover, method="pearson", na.rm=TRUE)
cor.test(transects$Evenness, transects$Richness, method="pearson", na.rm=TRUE)
cor.test(transects$Evenness, transects$Cover, method="pearson", na.rm=TRUE)
cor.test(transects$Cover, transects$Richness, method="pearson", na.rm=TRUE)

```

---

## Data transformations

```r
#Data transformations

#load libraries
library(car) #for QQ-plots
library(moments) #for analyzing skewness and kurtosis

#Let's look at the lab mouse temperature data to assess normality

#create numeric data
y <- as.numeric(c(97.9, 98.0, 98.2, 98.3, 98.3, 98.4, 98.4, 98.5, 98.5, 98.5,
                  98.6, 98.6, 98.6, 98.7, 98.7, 98.7, 98.8, 98.9, 98.9, 99.0,
                  99.1, 99.2, 99.2, 99.3, 99.4, 99.5, 99.7, 99.9))

#A histogram is one way to assess normality
hist(y, xlab = "Temperature", main = "") #looks pretty good
#What if we include more breaks?
hist(y, breaks = 8, xlab = "Temperature", main = "") #there's a bit of right skew
#We can also make a QQ plot
#The y-axis is the data's distribution and the x-axis is the expected normal distribution
#Look for a linear fit
qqnorm(y)
qqline(y) #most points are near the line
#Make a plot with the 95% confidence intervals
qqPlot(y, ylab = "Temperature") #most data are within the CLs 

#Analyze skewness and kurtosis 
skewness(y) #<1
kurtosis(y) #<3

#A Shapiro test is useful for assessing normality
#It tests the null hypothesis that the data is normally distributed
shapiro.test(y) #P = 0.698
#Seems like the data is normally distributed

#Look at the byssal threads data 
threads <- read.csv("Data/Byssal threads.csv")
hist(threads$Strength, xlab = "Strength", breaks = 8, main = "") #kind of right skewed
skewness(threads$Strength) #<1
kurtosis(threads$Strength) #<3
#Make a QQ-plot
qqPlot(threads$Strength, ylab = "Strength") #not great
#Try a Shapiro test
shapiro.test(threads$Strength) #P = 0.088
#It's safe to assume that normality is met, but can it be improved?

#Let's try a log transformation
threads$logstrength <- log(threads$Strength)
threads #notice logstrength is added to the data as a new column
hist(threads$logstrength, breaks = 8, xlab = "Log(Strength)",main = "") #kind of left skewed
shapiro.test(threads$logstrength) #P = 0.001
#That made it worse

#How about a square root transformation
threads$sqrtstrength <- sqrt(threads$Strength)
threads #sqrtstrength is another new column
hist(threads$sqrtstrength, xlab = "Sqrt(Strength)", breaks = 8, main = "") #not terrible
skewness(threads$sqrtstrength) #>-1
kurtosis(threads$sqrtstrength) #<3
qqPlot(threads$sqrtstrength, ylab = "Sqrt(Strength)")
shapiro.test(threads$sqrtstrength) #P = 0.22
#This is a slight improvement

#Look at kelp Bass gonad mass data
gonads <- read.csv("Data/Gonad data.csv")
hist(gonads$gonad_mass, xlab = "Gonad mass", breaks = 20, main = "") #a lot of right skew
skewness(gonads$gonad_mass) #>1
kurtosis(gonads$gonad_mass) #>3
qqPlot(gonads$gonad_mass, ylab = "Gonad mass")
shapiro.test((gonads$gonad_mass)) #P<<0.0001
#The data is far from normally distributed

#Try a log transformation
gonads$logmass <- log(gonads$gonad_mass)
hist(gonads$logmass, xlab = "Log(Gonad mass)", breaks = 20) #better but its left skewed
skewness(gonads$logmass) #>-1
kurtosis(gonads$logmass) #>3
qqPlot(gonads$logmass, ylab = "Log(Gonad mass)")
shapiro.test(gonads$logmass) #P<<0.0001
#A slight improvement

#Square root transformation
gonads$sqrtmass <- sqrt(gonads$gonad_mass)
hist(gonads$sqrtmass, xlab = "Sqrt(Gonad mass)", breaks = 20) #kind of better 
skewness(gonads$sqrtmass) #<1
kurtosis(gonads$sqrtmass) #<3
qqPlot(gonads$sqrtmass, ylab = "Sqrt(Gonad mass)")
shapiro.test(gonads$sqrtmass) #P<<0.0001
#Better than the raw data and the log transformation

```

---

## erodium snippets

```r
library(tidyverse)
library(car)
library(rstatix)

erodiumdata <- read.csv("Data/Erodium data.csv") %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("control", "shade", "wind", "3XH2O")
    ),
    plant_size = Width * Height,
    log_size = log(plant_size)
  )

# Check residual normality using the log-transformed model
log_model <- lm(log_size ~ Treatment, data = erodiumdata)
qqPlot(residuals(log_model))
shapiro.test(residuals(log_model))

# Test the treatment effect; Welch's test allows unequal variances
oneway.test(log_size ~ Treatment, data = erodiumdata, var.equal = FALSE)

# Compare each treatment pair; Gamesâ€“Howell allows unequal variances
games_howell_test(erodiumdata, log_size ~ Treatment)
```

---

## export_plot

```r
export_plot <- function(filename, plot = last_plot(), width = 7, height = 5) {
  if (!dir.exists("output")) dir.create("output")
  filepath <- file.path("output", filename)
  
  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(filepath, plot = plot, width = width, height = height, dpi = 300)
    message(paste("Saved ggplot to:", filepath))
  } else {
    dev.copy(jpeg, filename = filepath, width = width * 100, height = height * 100)
    dev.off()
    message(paste("Saved base R plot to:", filepath))
  }
}
```

---

## inspect_df

```r
inspect_df <- function(df, digits = 2) {
  # --- 1. Dataset Overview ---
  cat("========================================================\n")
  cat("                  DATASET OVERVIEW                      \n")
  cat("========================================================\n")
  cat(sprintf("Dimensions : %d rows x %d columns\n", nrow(df), ncol(df)))
  total_cells <- prod(dim(df))
  total_na <- sum(is.na(df))
  cat(sprintf("Total NA   : %d (%.1f%% of all cells)\n", total_na, (total_na / total_cells) * 100))
  cat(sprintf("Completeness: %d complete cases (%.1f%%)\n\n", 
              sum(complete.cases(df)), 
              (sum(complete.cases(df)) / nrow(df)) * 100))
  
  # Helper for statistical mode
  calc_mode <- function(x) {
    clean_x <- x[!is.na(x)]
    if (length(clean_x) == 0) return(NA)
    ux <- unique(clean_x)
    tab <- tabulate(match(clean_x, ux))
    ux[tab == max(tab)][1] # returns first mode if multimodal
  }
  
  # --- 2. Numeric Summary Table ---
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) > 0) {
    cat("--- Numeric Variables (Parametric & Non-Parametric) ---\n")
    num_list <- lapply(num_cols, function(col) {
      vals <- df[[col]]
      clean <- vals[!is.na(vals)]
      n_clean <- length(clean)
      
      if (n_clean > 1) {
        m    <- mean(clean)
        med  <- median(clean)
        mod  <- calc_mode(clean)
        s_d  <- sd(clean)
        vari <- var(clean)
        sem  <- s_d / sqrt(n_clean)
        iqr_val <- IQR(clean)
        mini <- min(clean)
        maxi <- max(clean)
      } else {
        m <- med <- mod <- s_d <- vari <- sem <- iqr_val <- mini <- maxi <- NA
      }
      
      data.frame(
        Variable = col,
        N = n_clean,
        NA_Count = sum(is.na(vals)),
        Mean = round(m, digits),
        SD = round(s_d, digits),
        SEM = round(sem, digits),
        Var = round(vari, digits),
        Median = round(med, digits),
        IQR = round(iqr_val, digits),
        Min = round(mini, digits),
        Max = round(maxi, digits),
        Mode = round(mod, digits),
        stringsAsFactors = FALSE
      )
    })
    print(do.call(rbind, num_list), row.names = FALSE)
    cat("\n")
  }
  
  # --- 3. Categorical / Factor Summary Table ---
  cat_cols <- names(df)[!sapply(df, is.numeric)]
  if (length(cat_cols) > 0) {
    cat("--- Categorical & Factor Variables ---\n")
    cat_list <- lapply(cat_cols, function(col) {
      vals <- df[[col]]
      clean <- vals[!is.na(vals)]
      top_val <- if (length(clean) > 0) names(sort(table(clean), decreasing = TRUE))[1] else NA
      top_freq <- if (length(clean) > 0) max(table(clean)) else 0
      
      data.frame(
        Variable = col,
        Type = class(vals)[1],
        Unique = length(unique(clean)),
        NA_Count = sum(is.na(vals)),
        Top_Category = as.character(top_val),
        Top_Freq = sprintf("%d (%.1f%%)", top_freq, (top_freq / max(1, length(clean))) * 100),
        stringsAsFactors = FALSE
      )
    })
    print(do.call(rbind, cat_list), row.names = FALSE)
  }
  cat("========================================================\n")
}



assumption_report <- function(data,
                              response = NULL,
                              group = NULL,
                              alpha = 0.05,
                              plots = TRUE) {
  
  # Confirm that the input is a data frame
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  
  # If no response column is specified, find all numeric columns
  if (is.null(response)) {
    response <- names(data)[vapply(data, is.numeric, logical(1))]
  }
  
  if (length(response) == 0) {
    stop("No numeric response columns were found.")
  }
  
  # Check that requested columns exist
  requested_columns <- c(response, group)
  
  if (!all(requested_columns %in% names(data))) {
    missing_columns <- requested_columns[!requested_columns %in% names(data)]
    stop(
      paste(
        "These columns were not found:",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  results <- list()
  
  for (y in response) {
    
    # Keep only the relevant columns
    if (is.null(group)) {
      working_data <- data[, y, drop = FALSE]
    } else {
      working_data <- data[, c(y, group), drop = FALSE]
      working_data[[group]] <- as.factor(working_data[[group]])
    }
    
    # Remove incomplete rows
    working_data <- working_data[complete.cases(working_data), ]
    values <- working_data[[y]]
    
    cat("\n====================================\n")
    cat("Response variable:", y, "\n")
    cat("Number of complete observations:", length(values), "\n")
    cat("====================================\n")
    
    if (length(values) < 3) {
      cat("Not enough observations for a Shapiro-Wilk test.\n")
      next
    }
    
    # Overall Shapiro-Wilk test
    if (length(values) <= 5000) {
      overall_shapiro <- shapiro.test(values)
      
      cat("\nOverall Shapiro-Wilk test:\n")
      print(overall_shapiro)
      
      if (overall_shapiro$p.value < alpha) {
        cat("Interpretation: evidence against normality.\n")
      } else {
        cat("Interpretation: no significant evidence against normality.\n")
      }
    } else {
      overall_shapiro <- NULL
      cat("\nShapiro-Wilk test skipped because n > 5000.\n")
    }
    
    # Results by group
    group_shapiro <- NULL
    bartlett_result <- NULL
    levene_result <- NULL
    
    if (!is.null(group)) {
      
      cat("\nShapiro-Wilk tests within each group:\n")
      
      group_values <- split(values, working_data[[group]])
      
      group_shapiro <- lapply(group_values, function(x) {
        if (length(x) >= 3 && length(x) <= 5000) {
          shapiro.test(x)
        } else {
          NULL
        }
      })
      
      for (group_name in names(group_shapiro)) {
        cat("\nGroup:", group_name, "\n")
        
        if (is.null(group_shapiro[[group_name]])) {
          cat("Not enough observations or n > 5000.\n")
        } else {
          print(group_shapiro[[group_name]])
        }
      }
      
      # Bartlett's test
      formula <- reformulate(group, response = y)
      bartlett_result <- bartlett.test(formula, data = working_data)
      
      cat("\nBartlett's test for equal variances:\n")
      print(bartlett_result)
      
      # Levene's test, if car is installed
      if (requireNamespace("car", quietly = TRUE)) {
        levene_result <- car::leveneTest(formula, data = working_data)
        
        cat("\nLevene's test for equal variances:\n")
        print(levene_result)
      } else {
        cat("\nLevene's test was skipped.\n")
        cat("Install the car package with: install.packages(\"car\")\n")
      }
    }
    
    # Q-Q plots
    if (plots) {
      if (is.null(group)) {
        qqnorm(
          values,
          main = paste("Q-Q Plot:", y)
        )
        qqline(values, col = "red", lwd = 2)
      } else {
        boxplot(
          values ~ working_data[[group]],
          main = paste("Boxplot:", y),
          xlab = group,
          ylab = y
        )
      }
    }
    
    # Save results
    results[[y]] <- list(
      response = y,
      n = length(values),
      overall_shapiro = overall_shapiro,
      group_shapiro = group_shapiro,
      bartlett = bartlett_result,
      levene = levene_result
    )
  }
  
  invisible(results)
}




#ANOVA-specific function
anova_assumption_report <- function(data,
                                    formula,
                                    alpha = 0.05,
                                    plots = TRUE) {
  
  # Fit the ANOVA model
  model <- lm(formula, data = data)
  
  # Extract residuals
  residual_values <- residuals(model)
  fitted_values <- fitted(model)
  
  # Remove missing or non-finite residuals
  residual_values <- residual_values[is.finite(residual_values)]
  
  cat("\n====================================\n")
  cat("ANOVA ASSUMPTION REPORT\n")
  cat("====================================\n")
  
  cat("\nModel:\n")
  print(formula)
  
  cat("\nSample size:", length(residual_values), "\n")
  
  # Shapiro-Wilk test on residuals
  if (length(residual_values) >= 3 &&
      length(residual_values) <= 5000) {
    
    shapiro_result <- shapiro.test(residual_values)
    
    cat("\nShapiro-Wilk test on model residuals:\n")
    print(shapiro_result)
    
    if (shapiro_result$p.value < alpha) {
      cat("Interpretation: evidence against normally distributed residuals.\n")
    } else {
      cat("Interpretation: no significant evidence against normal residuals.\n")
    }
    
  } else {
    shapiro_result <- NULL
    cat("\nShapiro-Wilk test could not be performed.\n")
  }
  
  # Bartlett and Levene tests
  bartlett_result <- NULL
  levene_result <- NULL
  
  if (requireNamespace("car", quietly = TRUE)) {
    
    cat("\nLevene's test for equal variances:\n")
    levene_result <- car::leveneTest(formula, data = data)
    print(levene_result)
    
  } else {
    
    cat("\nLevene's test was skipped.\n")
    cat("Install it with: install.packages(\"car\")\n")
  }
  
  # Diagnostic plots
  if (plots) {
    
    old_par <- par(mfrow = c(1, 2))
    on.exit(par(old_par), add = TRUE)
    
    # Q-Q plot of residuals
    qqnorm(
      residual_values,
      main = "Q-Q Plot of ANOVA Residuals"
    )
    qqline(residual_values, col = "red", lwd = 2)
    
    # Residuals versus fitted values
    plot(
      fitted_values,
      residuals(model),
      main = "Residuals vs. Fitted Values",
      xlab = "Fitted values",
      ylab = "Residuals"
    )
    abline(h = 0, col = "red", lwd = 2)
  }
  
  # Return results invisibly for later use
  invisible(list(
    model = model,
    residuals = residual_values,
    shapiro = shapiro_result,
    levene = levene_result
  ))
}
```

---

## log tryout with prof

```r
logcirc <- log(data$Circumference..cm.)
hist(logcirc, breaks = 30, main = "", xlab= "Tree height cm", xlim = c(0,10), ylim = c(0,25)) 
abline(v=lowerCL, col="blue")
abline(v=upperCL, col="blue")
```

---

## Logistic regression

```r
#Logistic regression

#load libraries
library(lme4)
library(car)
library(pscl) #for producing rho-squared
library(ggplot2)

#load data
perchdata <- read.csv("Data/Surfperch data.csv")
perchdata

#create generalized linear model
glm <- glm(mated ~ male.size, family = "binomial"(link="logit"), data=perchdata)
#look at the coefficients (estimate)
summary(glm)
#convert body mass to odds
exp(0.6563) #e^b1
#run log likelihood ratio test
Anova(glm)

#calculate rho-squared
pR2(glm)
#llh is the log likelihood from the fitted model
#llhnull is the log likelihood from the intercept-only model
#try calculating rho-square by hand
rho2 = 1 - (21.7156206/32.0517739)

#G2 method
#This does what Anova() does behind the scenes
G2 <- 2*((-21.7156206)-(-32.0517739)) #same as LR Chisq
1-pchisq(G2, df=1) #gives the P-value

#make simple plot
plot(perchdata$male.size, perchdata$mated, xlab="Male Size", ylab="Probability of Mating")
curve(predict(glm,data.frame(male.size=x),type="resp"),add=TRUE) # draws a curve

#make a nicer plot with ggplot()
predicted.data<-as.data.frame(predict(glm, type="response", se=TRUE))

ggplot(perchdata, aes(x=perchdata$male.size, y=perchdata$mated)) +
  theme_bw() +
  theme(axis.text.x=element_text(color="black", size=12), 
        axis.text.y=element_text(color="black", size=12), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
        geom_point() + geom_smooth(data = predicted.data, aes(y=fit)) +
        labs(x="Male Size (cm)", y="Pr(Mating Success)")

```

---

## Midterm Q2

```r
#Since the Q asks if density "differs", use a 2 tailed test 
#H_0: Average Buckwheat density is 4 plants per 25m^2 
#H_A: Average Buckwheat density differs from 4 plants per 25m^2
buckwheat <- read.csv("C:/Biostatistics/data/Buckwheat data.csv")

names(buckwheat)
str(buckwheat)

buckwheat_density <- buckwheat$Density

#Ttest does sample mean differ from 4
#t.test(buckwheat_density, mu = 4, alternative = "two.sided")
test <- t.test(
  buckwheat_density,
  mu = 4,
  conf.level = 0.95
)

test

#Boxplot of Buckwheat density
boxplot(
  buckwheat_density,
  notch = TRUE,
  ylab = "Eriogonum plants per 25m^2",
  main = "Buckwheat Density",
  col = "lightgreen"
)
#h adds horizontal H_0 line to boxplot
abline(
  h = 4,
  col = "red",
  lty = 2,
  lwd = 2
)

#p statistic < than 0.05, reject H_0. Buckwheat density is significantly greater than 4/25m^2
```

---

## Midterm Q3

```r
#Count data, so expected frequencies are fitted to a Poisson distribution 

#H_0: The data follows a Poisson distribution.
#Test of choice: Chi-Square Goodness of Fit Test
observed <- c(13, 26, 21, 13, 9)

expected <- c(12.69, 23.68, 22.09, 13.74, 9.79)

chisq.test(
  x = observed,
  p = expected / sum(expected),
  rescale.p = TRUE
)

#Poisson degrees of freedom = number of categories - 2, so 5 - 2 = 3
#because p-value is greater than 0.05 (~0.9831), we fail to reject the null hypothesis and conclude that the data follows a Poisson distribution.

```

---

## OrangeTree_Circumference_Script

```r
#Imports Spreadsheet
data <- read.csv("../Data/Orange trees FA26 - Sheet1.csv")

#Opens Spreadsheet In New Tab
View(data)

#Calculates Mean Circumference for Orange Trees
meancir <- data %>%
  summarise(mean(Circumference..cm.))

#Calculates Standard Deviation for Circumference of Orange Trees
stdcir <- data %>%
  summarise(sd(Circumference..cm.))

#Summary Function Above Generated Summary of Respective Data in a "mini-Table"

#Converts "mini-Table" Summary into Numeric Value
meancir <- as.numeric(meancir)

#Alt method: n <- sum(!is.na(data$Circumference..cm.))
#Alt method: SEMy <- sd(data$Circumference..cm., na.rm = TRUE) / sqrt(n)

#calculate confidence limits
#first calculate summary stats for the entire data set
ybar <- mean(data$Circumference..cm.)
SEMy <- sd(data$Circumference..cm.)/sqrt(count(data))

#Alt method: t_critical <- qt(0.975, df = n - 1)

#find the critical value of t with n-1 df
qt(p = 0.025, df = 131, lower.tail = FALSE)

lowerCL <- ybar-(SEMy*1.978239)
lowerCL
upperCL <- ybar+(SEMy*1.978239)
upperCL

#make a histogram w/ the confidence limits
hist(data$Circumference..cm., 
     breaks = 30, 
     main = "", 
     xlab= "Circumference", 
     ylim = c(0,25)
    ) 
abline(v=lowerCL, col="blue")
abline(v=upperCL, col="blue")

install.packages("car")
library(car)

qqnorm(data$Circumference..cm.) 
qqline(data$Circumference..cm., col = "red", lwd = 2)
boxplot(y, horizontal = TRUE, main = "Boxplot of response")

qqPlot(
  data$Circumference..cm.,
  distribution = "norm",
  envelope = 0.95,
  id = FALSE,
  main = "Normal Q-Q Plot of Tree Circumference"
)

```

---

## Paired t-test

```r
#Paired t-test

#load libraries
library(tidyverse)
library(moments)
library(ggpubr)

#perform paired t-test on lab mice data

#create numeric data
#before drug
x <- as.numeric(c(98.5, 98.7, 98.8, 98.7, 98.6, 98.7, 98.7, 98.4, 98.5, 98.7,
                  98.2, 98.3, 98.5, 98.9, 98.2, 98.9, 98.8, 98.6, 98.4, 98.9, 
                  98.5, 98.5, 98.9, 98.4, 98.4, 98.7, 98.4, 98.7)) 
#after drug
y <- as.numeric(c(97.9, 98.0, 98.2, 98.3, 98.3, 98.4, 98.4, 98.5, 98.5, 98.5,
                  98.6, 98.6, 98.6, 98.7, 98.7, 98.7, 98.8, 98.9, 98.9, 99.0,
                  99.1, 99.2, 99.2, 99.3, 99.4, 99.5, 99.7, 99.9))

#combine samples into a data frame
drugdata <- data.frame(Before = x, After = y)
drugdata

#calculate differences
drugdata$Differences <- drugdata$Before - drugdata$After
drugdata #the above added a new column to the data frame

#calculate summary stats
mean <- mean(drugdata$Differences)
mean
stdev <- sd(drugdata$Differences)
stdev
n <- count(drugdata)
n
SEM <- stdev/sqrt(n) #SEM with n=28
SEM

#check model assumption of normality
hist(drugdata$Differences, breaks = 6, main = '') #looks good
hist(drugdata$Differences, breaks = 8, main = '') #bit of a negative skew
qqPlot(drugdata$Differences) #not too bad
skewness(drugdata$Differences) # >-1 and <1
kurtosis(drugdata$Differences) # <3
shapiro.test(drugdata$Differences) #P > 0.05

#calculate t
t <- (mean-0)/SEM
t

#perform test
#one tailed test (only concerned with temperatures elevating)
t.test(drugdata$Differences, alternative = c("less"), mu = 0, paired = FALSE)
#two-tailed test (temperatures could have gone in either direction)
t.test(drugdata$Differences, mu = 0, paired = FALSE)


#paired t-test on for chaparral abundance data

#load data
chaparraldata <- read.csv("Data/Chaparral data.csv")
chaparraldata
#to load data into a new window 
view(chaparraldata)

#check normality
hist(chaparraldata$Differences, main = '') #looks normal-ish
qqPlot(chaparraldata$Differences) #looks ok
skewness(chaparraldata$Differences) # <1
kurtosis(chaparraldata$Differences) # <3
shapiro.test(chaparraldata$Differences) #P > 0.05

#calculate summary stats
mean <- mean(chaparraldata$Differences)
mean
stdev <- sd(chaparraldata$Differences)
stdev
n <- count(chaparraldata)
n
SEM <- stdev/sqrt(n) #n=9
SEM

#calculate t
t <- (mean-0)/SEM #hypothesized mean difference of pairs = 0
t

#perform t-test
#two tailed test
t.test(chaparraldata$Desert, chaparraldata$Coast, paired = TRUE)
#one-tailed test (only tests if the abundances are higher near the coast)
t.test(chaparraldata$Desert,  alternative = c("less"), chaparraldata$Coast, paired = TRUE)

#create a paired box plot

#transform into long data
#the %>% operator is called a pipe - it tells R to take the data and then gather it
chaplong <- chaparraldata %>%
  gather(key = "site", value = "abundance", Desert, Coast)
chaplong #notice how the data is arranged now

#use ggpubr for plotting
ggpaired(chaplong, x = "site", y = "abundance", 
                order = c("Desert", "Coast"),
                ylab = "Relative abundance", xlab = "")

#or use boxplot() without lines
boxplot(chaparraldata$Coast, chaparraldata$Desert, names = c("Coast", "Desert"),
                ylab="Relative abundance")

```

---

## problem 4

```r
#2 independent groups, Wildtype and Mutant, are being compared for their chlorophyll content. Rows arent matched pairs, so an independent 2 sample comparison
#H_0: The means of Chlorophyll in the two groups are equal
#H_A: The means of Chlorophyll in the two groups are not equal
#First observe the data and calculate the means and sums, sample size, & standard error of the two groups. 
#Then, we will check the assumptions of normality and equal variance before deciding on the appropriate statistical test to use.

library(tidyverse) #for data wrangling
library(ggplot2) #for plotting
library(car) #for qqPlot function
library(moments)

Microalgaedata <- read.csv("data/Microalagae data.csv")
Microalgaedata
inspect_df(Microalgaedata)

Microalgaedata_long <- Microalgaedata %>%
  gather(key = "Genotype", value = "Concentration", Wildtype, Mutant, na.rm = TRUE)
Microalgaedata_long
inspect_df(Microalgaedata_long)

qqPlot(Microalgaedata_long$Concentration, ylab = "Concentration", distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Concentration Q-Q plot", col.lines = "red")
shapiro.test(Microalgaedata_long$Concentration) #H_0: The data is normally distributed. H_A: The data is not normally distributed. If p-value is less than 0.05, the data is not normally distributed
bartlett.test(Concentration ~ Genotype, data = Microalgaedata_long)

Microalgaedata_long$log_Concentration <- log(Microalgaedata_long$Concentration)
Microalgaedata_long
qqPlot(Microalgaedata_long$log_Concentration, ylab = "Log Concentration", distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Log Concentration Q-Q plot", col.lines = "red")
shapiro.test(Microalgaedata_long$log_Concentration)  #H_0: The data is normally distributed. H_A: The data is not normally distributed. If p-value is less than 0.05, the data is not normally distributed
bartlett.test(log_Concentration ~ Genotype, data = Microalgaedata_long)

hist(Microalgaedata_long$Concentration, main = "Histogram of Concentration", xlab = "Concentration", col = "lightblue", border = "black")
hist(Microalgaedata_long$log_Concentration, main = "Histogram of Log Concentration", xlab = "Log Concentration", col = "lightblue", border = "black")
boxplot(Microalgaedata$Wildtype, Microalgaedata$Mutant, names = c("Wildtype", "Mutant"),
        col = "skyblue", xlab = "Strain", ylab = "Concentration")
boxplot(log_Concentration ~ Genotype, data = Microalgaedata_long, col = "skyblue", xlab = "Strain", ylab = "Log Concentration")
t.test(log_Concentration ~ Genotype, data = Microalgaedata_long, var.equal = TRUE) 

```

---

## PSet2Q1Resubmit

```r
#H_0: The sex ratio of Hymenoptera Bees is equal to 1:1
#H_A: The workerbees are tipping the sex ratio to 3:1

#Goodness-of-fit tests

#Chi-square test for equal  sex ratio in Hymenoptera Bees

#create data
female<-371
male<-250

#set up table of observed and expected value:
observed <- matrix(c(female, male), nrow=1, ncol=2) #created a table of observed values
exp <- matrix(c((female+male)*0.5, (female+male)*0.5), nrow=1, ncol=2)

#create a matrix to view data 
beedata <- matrix(c(observed, exp), nrow=2, ncol=2)
beedata
#add row and column names
row.names(beedata) <- (c("Female", "Male"))
colnames(beedata) <- (c("Observed", "Expected"))
beedata

#perform X2 test:
X2test<- sum((observed-exp)^2/exp)
X2test #chi-square test statistic, df = #categories-1=1
1-pchisq(X2test, df=1) #P

#perform G2 test
#calculate G:
G2 <- 2*sum(observed*log(observed/exp))
G2 
#get p-value
1-pchisq(G2, df=1) #P 

#X2 test and G2 test are both significant(~24, ~24 respectively), so we reject the null hypothesis that the sex ratio of Hymenoptera Bees is equal to 1:1. 
#The data suggests that the worker bees are tipping the sex ratio to 3:1.

```

---

## PSet2Q2Resubmit

```r
#H_0: The observed data fits the expected 9:3:3:1 ratio
#H_A: The observed data does not fit the expected 9:3:3:1 ratio
#Critical value: 7.815
data <- read.csv("Data/Dihybrid Cross Gtest - Sheet1.csv")
view(data)
#Goodness-of-fit tests

#create Observed data
yellowround<-96
yellowwrinkled<-18
greenround<-24
greenwrinkled<-7

n<-sum(yellowround, yellowwrinkled, greenround, greenwrinkled)  #created a table of observed values
exp<-matrix(c((n*9/16), (n*3/16), (n*3/16), (n*1/16)), nrow=1, ncol=4)
observed <- matrix(c(yellowround, yellowwrinkled, greenround, greenwrinkled), nrow=1, ncol=4) 
#created a table of observed values

#create a matrix to view data 
peadata <- matrix(c(observed, exp), nrow=4, ncol=2)
peadata
#add row and column names
row.names(peadata) <- (c("yellowround", "yellowwrinkled", "greenround", "greenwrinkled"))
colnames(peadata) <- (c("Observed", "Expected"))
peadata            
            
#set up table of observed and expected value:

#perform G2 test
#calculate G:
G2 <- 2*sum(observed*log(observed/exp))
G2 #Observe G value is less than critical 7.815
#get p-value
p_value <- pchisq(G2, df = 3, lower.tail = FALSE)
p_value

#Conclusion: Since the p-value is greater than 0.05 (0.077), we fail to reject the null hypothesis 
#and conclude that there's not enough evidence say they differ from expected 9:3:3:1 ratio. 
```

---

## PSet2Q4Resubmit

```r
#H_0: The observed data fits the a priori expected 95% of trees are expected to produce fruit
#H_A: The observed data does not fit the a priori expected 95% of trees are expected to produce fruit
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

#Conclusion: since all 132 trees had fruit, the one sided binomial test did not 
#show evidence that the fruiting tree proportion is less than 95% (p=1!)
```

---

## PSet2Q5Resubmit

```r
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
```

---

## PSet3Q2

```r
load_biostats_packages()
#Contingency tables
#Imports Spreadsheet
data <- read.csv("Data/Ladybird data.csv")

#Opens Spreadsheet In New Tab
View(data)

#load libraries 
library(tidyverse)
library(DescTools) #for G test
library(ggplot2) #for plotting

#analyze dispersal and reproduction data

#create data
redrural <- 70
redindustrial <- 85
blackrural <- 30
blackindustrial <- 115

# Make contingency table using ONLY the observed values
contable_obs <- matrix(c(redrural, blackrural,
                         redindustrial, blackindustrial),
                       nrow = 2, byrow = TRUE)

# Name rows and columns
rownames(contable_obs) <- c("Rural", "Industrial")
colnames(contable_obs) <- c("Red", "Black")

addmargins(contable_obs)

# Calculate expected values
contable_exp <- outer(rowSums(contable_obs),
                      colSums(contable_obs)) / sum(contable_obs)

contable_exp

# Perform G-test
GTest(contable_obs)

#make plot
names <- matrix(c("Rural", "Ind"))
barplot(contable_obs, beside = TRUE, legend.text = rownames(contable_obs), xlab = "Color", ylab = "Counts", args.legend = list(x = "topright"))

```

---

## PSet3Q3 residual snippets not final

```r
#Create Data
ConeFlying <- 
  BarkFlying <-
  LeavesFlying <-
  BranchFlying  <-
  
  flyingcone <- data %>%
  filter(Behavior == "Flying")

Countdata <- function(data, x, y){
  countdata <- data %>%
    filter(behavior == x) %>%
    filter(Substrate == y) %>%
    summarise(Count == sum(Count))
  unlist(Countdata)
}

make_contingency_table <- function(data) {
  xtabs(Count ~ Substrate + Behavior, data = data)
}

contable_obs <- make_contingency_table(data)
contable_obs

contable_obs <- matrix(c())

library(ggplot2)

ggplot(data, aes(x = Substrate, y = Count, fill = Behavior)) + 
  geom_col(position = "dodge") + 
  labs(
    title = "Foraging Behaviors of Mountain Chickadees by Substrate",
    x = "Substrate",
    y = "Observed Count",
    fill = "Behavior"
  ) + 
  theme_minimal()





fisher.test(contable_obs, simulate.p.value=TRUE)

library(ggplot2)


ConeFlying
2	Bark	Flying
1	Leaves	Flying
1	Branch	Flying
0	Cone	Hopping
0	Bark	Hopping
13	Leaves	Hopping
23	Branch	Hopping
0	Cone	Climbing
1	Bark	Climbing
0	Leaves	Climbing
0	Branch	Climbing
3	Cone	Gleaning
0	Bark	Gleaning
10	Leaves	Gleaning
5	Branch	Gleaning
2	Cone	Probing
0	Bark	Probing
0	Leaves	Probing
10	Branch	Probing
ConeHanging
BarkHanging
LeavesHanging
BranchHanging
```

---

## PSet3Q3

```r
#load libraries 
load_biostats_packages()
library(tidyverse)
library(DescTools) #for G test
library(ggplot2) #for plotting

# CSV import
data <- read.csv("data/Foraging data.csv")
view(data)

#make_contingency_table <- function_ Creates a function named make_contingency_table. Data; dataset containing my observations or summarized counts. count_col; column containing the numerical counts, such as "Count". row_col; categorical variable I want as the table rows, such as "Substrate".col_colThe categorical variable you want as the table columns, such as Behavior.
make_contingency_table <- function(data, count_col, row_col, col_col) {
#This section builds the formula; table_formula <- reformulate(c(row_col, col_col),response = count_col) 
#reformulate() builds an R formula using the column names you provide  
    table_formula <- reformulate(
    c(row_col, col_col),
    response = count_col
  )
#xtabs is an R function that creates a cross-tabulation, aka a contingency table
#The first data tells xtabs which argument it is receiving, The second data refers to our actual dataset in the Global Environment.
        xtabs(table_formula, data = data)
}

# Apply it to the current dataset
#Runs your make_contingency_table function, Uses the dataset named data, Saves the resulting contingency table as contable_obs
contable_obs <- make_contingency_table(
  data = data,
  count_col = "Count",
  row_col = "Substrate",
  col_col = "Behavior"
)

contable_obs

#run the appropriate analysis. [G-test?] 
#can't use G-test for this problem since we have zero values 
#DescTools::GTest(contable_obs)

#Non-parametric fisher test of independence template
#defining the function, B sets number of simulations
run_fisher_test <- function(conttingency_table, B = 1000) {
#runs fisher's exact test
    fisher.test(
    conttingency_table, 
#tells r to estimate p-value using random simulations instead of calculating it directly
        simulate.p.value = TRUE, 
    B = B
  )
}

#using it 
fisher_result <- run_fisher_test(contable_obs)

fisher_result

#creating template graph
plot_contingency <- function(
    data,
    count_col,
    x_col,
    fill_col,
    graph_type = "counts",
    graph_title = NULL
) {
  
# Allowing only these two graph types
  graph_type <- match.arg(graph_type, c("counts", "proportion"))
  
# Creates default title if not provided
  if (is.null(graph_title)) {
    graph_title <- paste(fill_col, "by", x_col)
  }
  
#Choose appropriate bar arrangement
#graph_type=counts; creates raw-count graph  
  bar_position <- if (graph_type == "counts") { 
#position=dodge creates side by side bars, best for comparing observed counts    
    "dodge"
  } else { 
#position=fill creates stacked bars, where every substrate totals 100%    
    "fill"
  }
  
#Creates graph
  graph <- ggplot(
    data,
#maps supplied collumns to graph. x_col; goes on x-axis. count_col; determines bar height. fill_col; determines bar
        aes(
#.data[[...]] lets function use column names supplied as text          
      x = .data[[x_col]],
      y = .data[[count_col]],
      fill = .data[[fill_col]]
    )
  ) + 
#geom_col; creates bars using existing values in data count column  
    geom_col(position = bar_position) +
    labs(
      title = graph_title,
      x = x_col,
      y = if (graph_type == "counts") {
        "Observed Count"
      } else { 
#graph_type=proportion creates relative frequency graph        
        "Proportion"
      },
      fill = fill_col
    ) +
    theme_minimal()
  
#Formats y-axis as percentages for proportional graphs
  if (graph_type == "proportion") {
    graph <- graph +
      scale_y_continuous(labels = scales::percent)
  }
  
  return(graph)
}

#using it with current data
#for future contingency problems, change only column name arguments: data=new_data, count_col=Frequency, x_col=Habitat, fill_col=Species, graph_type=counts
plot_contingency(
  data = data,
  count_col = "Count",
  x_col = "Substrate",
  fill_col = "Behavior",
  graph_type = "counts", #for proportional graph change to graph_type = "proportion"
  graph_title = "Mountain Chickadee Foraging Behaviors by Substrate"
)


```

---

## PSet3Q4 - w problem type template at end

```r
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
inspect_df(nectar)

#create histogram 
par(mfrow = c(1, 2))

hist(
  nectar,
  main = "Original Nectar Volumes",
  xlab = "Nectar volume (ÂµL)",
  col = "lightblue",
  border = "black"
)

boxplot(
  nectar,
  main = "Original Nectar Volumes",
  ylab = "Nectar volume (ÂµL)",
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
  xlab = "Nectar volume (ÂµL)",
  col = "lightblue"
)

boxplot(
  nectar,
  main = "Original Boxplot",
  ylab = "Nectar volume (ÂµL)",
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
```

---

## PSet3Q4

```r
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
  xlab = "Nectar volume (ÂµL)",
  col = "lightblue",
  border = "black"
)

boxplot(
  nectar,
  main = "Original Nectar Volumes",
  ylab = "Nectar volume (ÂµL)",
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
  xlab = "Nectar volume (ÂµL)",
  col = "lightblue"
)

boxplot(
  nectar,
  main = "Original Boxplot",
  ylab = "Nectar volume (ÂµL)",
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
```

---

## PSet3Q5

```r
# Import CSV file
protein_data <- read.csv("data/Protein data.csv")
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

qqPlot(protein_data$Concentration, ylab = "Protein Concentration (mg protein/mg tissue)", main = "Q-Q Plot of Protein Concentration")
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
view(protein_data)

```

---

## PSet4 Q3

```r
#H_0: The Abundance of Pacific Purple Sea Urchins is the same in both Rocky & Flat Substrates
#H_A: The Abundance of Pacific Purple Sea Urchins Differ Between Rocky & Flat Substrates
  
Urchins <- read.csv("data/PSet4 Q3 - Sheet1.csv")
view(Urchins)
inspect_df(Urchins)
VarianceFlat  <- var(Urchins$Flat)
VarianceRocky <- var(Urchins$Rocky)

CoeficientofVariationFlat <- sqrt(VarianceFlat) / mean(Urchins$Flat)
CoeficientofVariationRocky <- sqrt(VarianceRocky) / mean(Urchins$Rocky)

Ttest <- t.test(Urchins$Flat, Urchins$Rocky, alternative = "two.sided")
Ttest

#Conclusion: Our T-test statistic returned a P-value of 0.05622 which is just above the threshhold
#of 0.05, suggesting marginally significant difference in Pacific Purple Sea Urchin between 
#the Rocky and Flat Substrates
```

---

## PSet4 Q4

```r
#2 independent groups, Wildtype and Mutant, are being compared for their chlorophyll content. Rows arent matched pairs, so an independent 2 sample comparison
#H_0: The means of Chlorophyll in the two groups are equal
#H_A: The means of Chlorophyll in the two groups are not equal
#First observe the data and calculate the means and sums, sample size, & standard error of the two groups. 
#Then, we will check the assumptions of normality and equal variance before deciding on the appropriate statistical test to use.

library(tidyverse) #for data wrangling
library(ggplot2) #for plotting
library(car) #for qqPlot function

Microalgaedata <- read.csv("data/Microalagae data.csv")
Microalgaedata

mean(Microalgaedata$Wildtype)
mean(Microalgaedata$Mutant, na.rm = TRUE) #na.rm mean remove missing values
sum(Microalgaedata$Wildtype)
sum(Microalgaedata$Mutant, na.rm = TRUE)

sample.size <- c(length(Microalgaedata$Wildtype), length(na.omit(Microalgaedata$Mutant)))
sample.size

standard.error <- c(sd(Microalgaedata$Wildtype)/sqrt(length(Microalgaedata$Wildtype)), sd(na.omit(Microalgaedata$Mutant))/sqrt(length(na.omit(Microalgaedata$Mutant))))
standard.error

qqPlot(Microalgaedata$Wildtype, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Wildtype Q-Q plot", col.lines = "red")
shapiro.test(Microalgaedata$Wildtype) #H_0: The data is normally distributed. H_A: The data is not normally distributed. If p-value is less than 0.05, the data is not normally distributed

qqPlot(na.omit(Microalgaedata$Mutant), distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Mutant Q-Q plot", col.lines = "red")
shapiro.test(na.omit(Microalgaedata$Mutant))

                      #data.frame creates table with 2 columns, strain {a label for each measurement} and chlorophyll {the actual measurement values}
Microalgaedata_long <- data.frame(
  strain = c(rep("Wildtype", length(Microalgaedata$Wildtype)), rep("Mutant", length(na.omit(Microalgaedata$Mutant)))),
  chlorophyll = c(Microalgaedata$Wildtype, na.omit(Microalgaedata$Mutant)))
Microalgaedata_long

#Bartlett test H0: The variances are equal
bartlett.test(chlorophyll ~ strain, data = Microalgaedata_long)

Microalgaedata_long$log_chlorophyll <- log(Microalgaedata_long$chlorophyll)
qqPlot(log_chlorophyll ~ strain, data = Microalgaedata_long, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Log-transformed Q-Q plot", col.lines = "red") # "==" means "is equal to?" true or false answers. == asks, <- assigns
shapiro.test(Microalgaedata_long$log_chlorophyl[Microalgaedata_long$strain == "Wildtype"]) #Square brackets tell R which entries to take out of the log_chlorophyll column to include only those in the shapiro test. In this case, we are taking out the log_chlorophyll values for the Wildtype strain.
shapiro.test(Microalgaedata_long$log_chlorophyll[Microalgaedata_long$strain == "Mutant"]) #
bartlett.test(log_chlorophyll ~ strain, data = Microalgaedata_long) #Bartlett test H0: The variances are equal
#No evidence the variances differ!

#H_0: The mean Log of Chlorophyll concentration in the two groups are equal
#H_A: The mean Log of Chlorophyll concentration in the two groups are not equal
t.test(log_chlorophyll ~ strain, data = Microalgaedata_long, var.equal = TRUE) 
#df = n_Wildtype + n_Mutant - 2 = 30 + 25 - 2 = 53
#p-value = 0.1193 > 0.05, so we fail to reject the null hypothesis.
#There is not enough evidence to conclude that the mean Log of Chlorophyll concentration is different between the two groups.
```

---

## PSet4 Q5

```r
#H_0: individual Yew conifers produce roughly equal proportions of male and female cones (ie. they're monoecious), population mean index is 0.5
#H_A: individual Yew conifers do not produce equal proportions of male and female cones (ie. they're dioecious), population mean index is not 0.5
#I wanna do a Paired T-test because I'm looking at Dependent "paired" data, the Male and Female cone counts are from the same individual trees.

library(tidyverse)
library(ggplot2)
library(car)
library(readxl)

yewsdata <- read.csv("data/Yew data.csv")
inspect_df(yewsdata)
yewsdata

yewsdata_long <- yewsdata %>%
  gather(key = "Gender", value = "Abundance", Male, Female, na.rm = TRUE)
yewsdata_long

qqPlot(yewsdata_long$Abundance, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Abundance of Cones Q-Q plot", col.lines = "red")  

yewsdata_long$Abundance <- sqrt(yewsdata_long$Abundance)
yewsdata_long
qqPlot(sqrt(yewsdata_long$Abundance), distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Sqrt Abundance of Cones", col.lines = "red")

#yewslog_Male <- log(yewsdata$Male)
#qqPlot(yewslog_Male, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Male Cones LOG Q-Q plot", col.lines = "red") #can't graph Infinite values

wilcox.test(yewsdata$Male, yewsdata$Female, paired = TRUE, alternative = "two.sided", conf.level = 0.95)       

hist(yewsdata$PhGenRto, main = "Histogram of Phenotypic Gender Expression", xlab = "Phenotypic Gender Ratio", ylim = c(0, 35), ylab = "Frequency", col = "lightblue", border = "black")


```

---

## PSet4Q1.1 (2)

```r
#Imports Spreadsheet
orangedata <- read.csv("../Data/Orange trees FA26 - Sheet1.csv")

#Opens Spreadsheet In New Tab
View(orangedata)
orangedata
#orangedata2 <- orangedata[,-c(1,2,4,5,6,7)]
#orangedata2

#Calculates Mean Circumference for Orange Trees
meancir <- orangedata %>%
  summarise(mean(Circumference..cm.))
#Alt cleaner version: meancir <- mean(data$Circumference..cm., na.rm = TRUE) (na.rm = TRUE ignores missing values)
#Using Alt cleaner version negates the need to convert "mini-Table" Summary into Numeric Value below

#Calculates Standard Deviation for Circumference of Orange Trees
stdcir <- orangedata %>%
  summarise(sd(Circumference..cm.))
#Alt cleaner version: stdcir <- sd(orangedata$Circumference..cm., na.rm = TRUE)

#Summary Function Above Generated Summary of Respective Data in a "mini-Table"

#Converts "mini-Table" Summary into Numeric Value
meancir <- as.numeric(meancir)

#calculate confidence limits
#first calculate summary stats for the entire data set
ybar <- mean(orangedata$Circumference..cm.) #sample mean

#n <- sum(!is.na(orangedata$Circumference..cm.))
#SEMy <- sd(orangedata$Circumference..cm., na.rm = TRUE) / sqrt(n)
SEMy <- sd(orangedata$Circumference..cm.)/sqrt(count(orangedata)) #Standard Error of the Mean: SE = s/sqrt(n), where s = sample standard deviation and n = sample size

#find the critical value of t with n-1 df
#t_critical <- qt(0.975, df = n - 1) : better & reusable way to find critical t value for 95% confidence interval
qt(p = 0.025, df = 131, lower.tail = FALSE)
lowerCL <- ybar-(SEMy*1.978239) #This calculates:[lower limit = ybar - t*SE]
lowerCL
upperCL <- ybar+(SEMy*1.978239) #This calculates:[upper limit = ybar + t*SE]
upperCL

#complete reusable code for calculating confidence limits:
#{values <- data$Circumference..cm.
#values <- values[!is.na(values)]
#n <- length(values)
#ybar <- mean(values)
#s <- sd(values)
#SE <- s / sqrt(n)
#t_critical <- qt(0.975, df = n - 1)
#lowerCL <- ybar - t_critical * SE
#upperCL <- ybar + t_critical * SE
#lowerCL
#upperCL}

#make a histogram w/ the confidence limits
#hist(
#  orangedata$Circumference..cm., 
#  breaks = 30, 
#  main = "", 
#  xlab= "Circumference", 
#  ylim = c(0,25)
#  ) 
#abline(v=lowerCL, col="blue") #adds vertical line to histogram at lower confidence limit
#abline(v=upperCL, col="blue") #adds vertical line to histogram at upper confidence limit

#Compare quantiles of the data to the quantiles of a normal distribution
#qqnorm(orangedata$Circumference..cm.)
#qqline(orangedata$Circumference..cm., col = "red", lwd = 2)

# Select circumference values
circumference <- orangedata$Circumference..cm.

# Remove missing values
circumference <- circumference[!is.na(circumference)]

# Create normal Q-Q plot
qqnorm(
  circumference,
  main = "Normal Q-Q Plot of Tree Circumference",
  xlab = "Theoretical Normal Quantiles",
  ylab = "Observed Circumference"
)

# Add reference line
qqline(
  circumference,
  col = "red",
  lwd = 2
)

qqPlot(
  circumference,
  distribution = "norm", #Compares the quantiles of the data to the quantiles of a normal distribution
  envelope = 0.95, #Adds a 95% confidence envelope to the plot
  id = FALSE, #Prevents R from labeling outliers on the plot
  main = "Normal Q-Q Plot of Tree Circumference" #Adds a title to the plot
)
```

---

## PSet4Q1.1

```r
#Imports Spreadsheet
orangedata <- read.csv("../Data/Orange trees FA26 - Sheet1.csv")

#Opens Spreadsheet In New Tab
View(orangedata)
orangedata
#orangedata2 <- orangedata[,-c(1,2,4,5,6,7)])

#orangedata2

#Calculates Mean Circumference for Orange Trees
meancir <- orangedata %>%
  summarise(mean(Circumference..cm.))
#Alt cleaner version: meancir <- mean(data$Circumference..cm., na.rm = TRUE) (na.rm = TRUE ignores missing values)
#Using Alt cleaner version negates the need to convert "mini-Table" Summary into Numeric Value below

#Calculates Standard Deviation for Circumference of Orange Trees
stdcir <- orangedata %>%
  summarise(sd(Circumference..cm.))
#Alt cleaner version: stdcir <- sd(orangedata$Circumference..cm., na.rm = TRUE)

#Summary Function Above Generated Summary of Respective Data in a "mini-Table"

#Converts "mini-Table" Summary into Numeric Value
meancir <- as.numeric(meancir)

#calculate confidence limits
#first calculate summary stats for the entire data set
ybar <- mean(orangedata$Circumference..cm.) #sample mean

#n <- sum(!is.na(orangedata$Circumference..cm.))
#SEMy <- sd(orangedata$Circumference..cm., na.rm = TRUE) / sqrt(n)
SEMy <- sd(orangedata$Circumference..cm.)/sqrt(count(orangedata)) #Standard Error of the Mean: SE = s/sqrt(n), where s = sample standard deviation and n = sample size

#find the critical value of t with n-1 df
#t_critical <- qt(0.975, df = n - 1) : better & reusable way to find critical t value for 95% confidence interval
qt(p = 0.025, df = 131, lower.tail = FALSE)
lowerCL <- ybar-(SEMy*1.978239) #This calculates:[lower limit = ybar - t*SE]
lowerCL
upperCL <- ybar+(SEMy*1.978239) #This calculates:[upper limit = ybar + t*SE]
upperCL

#complete reusable code for calculating confidence limits:
#{values <- data$Circumference..cm.
#values <- values[!is.na(values)]
#n <- length(values)
#ybar <- mean(values)
#s <- sd(values)
#SE <- s / sqrt(n)
#t_critical <- qt(0.975, df = n - 1)
#lowerCL <- ybar - t_critical * SE
#upperCL <- ybar + t_critical * SE
#lowerCL
#upperCL}

#make a histogram w/ the confidence limits
#hist(
#  orangedata$Circumference..cm., 
#  breaks = 30, 
#  main = "", 
#  xlab= "Circumference", 
#  ylim = c(0,25)
#  ) 
#abline(v=lowerCL, col="blue") #adds vertical line to histogram at lower confidence limit
#abline(v=upperCL, col="blue") #adds vertical line to histogram at upper confidence limit

#Compare quantiles of the data to the quantiles of a normal distribution
#qqnorm(orangedata$Circumference..cm.)
#qqline(orangedata$Circumference..cm., col = "red", lwd = 2)

# Select circumference values
circumference <- orangedata$Circumference..cm.

# Remove missing values
circumference <- circumference[!is.na(circumference)]

# Create normal Q-Q plot
qqnorm(
  circumference,
  main = "Normal Q-Q Plot of Tree Circumference",
  xlab = "Theoretical Normal Quantiles",
  ylab = "Observed Circumference"
)

# Add reference line
qqline(
  circumference,
  col = "red",
  lwd = 2
)

qqPlot(
  circumference,
  distribution = "norm", #Compares the quantiles of the data to the quantiles of a normal distribution
  envelope = 0.95, #Adds a 95% confidence envelope to the plot
  id = FALSE, #Prevents R from labeling outliers on the plot
  main = "Normal Q-Q Plot of Tree Circumference" #Adds a title to the plot
)
```

---

## PSet4Q1

```r
# Load packages
library(car)
library(moments)

#Imports Spreadsheet
orangedata <- read.csv("Data/Orange trees FA26 - Sheet1.csv")

#Opens Spreadsheet In New Tab
View(orangedata)
orangedata #isolate collumns of interest: Circumference and Height of Trees
orangedata2 <- orangedata[,-c(1,2,4,5,6,7)]
orangedata2

#Calculate & renames Mean then Standard Deviation of Circumference for Orange Trees
meancir <- mean(orangedata2$Circumference..cm.)
meansd <- sd(orangedata2$Circumference..cm.)

n <- length(orangedata2$Circumference..cm.)
SEMmeancir <- meansd/sqrt(n) #Standard Error of the Mean: SE = s/sqrt(n), where s = sample standard deviation and n = sample size

#find the critical value of t with n-1 df
t_critical <- qt(1-0.05/2, df = n - 1) #: better & reusable way to find critical t value for 95% confidence interval
lowerCL <- meancir-(t_critical*SEMmeancir) #This calculates:[lower limit = ybar - t*SE]
lowerCL
upperCL <- meancir+(t_critical*SEMmeancir) #This calculates:[upper limit = ybar + t*SE]
upperCL

#Compare quantiles of the data to the quantiles of a normal distribution
qqPlot(orangedata2$Circumference, distribution = "norm", envelope = 0.95, main = "Normal Q-Q Plot of Tree Circumference")
skewness(orangedata2$Circumference..cm.) #near 0 means more symmetry, positive means longer right tail
kurtosis(orangedata2$Circumference..cm.) #normal distribution has kurtosis of 3 
shapiro.test(orangedata2$Circumference..cm.) #Shapiro p value greater than 0.5: no strong evidence against normality 

log_circumnference <- log(orangedata2$Circumference..cm.) #log transformation of circumference data
qqPlot(log_circumnference, distribution = "norm", envelope = 0.95, main = "Q-Q Plot of Log Circumference")
skewness(log_circumnference)
kurtosis(log_circumnference)
shapiro.test(log_circumnference)

sqrt_circumnference <- sqrt(orangedata2$Circumference..cm.) #square root transformation of circumference data
qqPlot(sqrt_circumnference, distribution = "norm", envelope = 0.95, main = "Q-Q Plot of Square Root Circumference")
skewness(sqrt_circumnference)
kurtosis(sqrt_circumnference)
shapiro.test(sqrt_circumnference)
#---------------------------------------------------------------
#Now for Height of Trees
height <- orangedata2$Tree.height..cm.
height <- height[!is.na(height)]

meanheight <- mean(height)
sdheight <- sd(height)
nheight <- sum(!is.na(height))

SEMheight <- sdheight/sqrt(nheight)
t_critical_height <- qt(0.975, df = nheight - 1)

lower_height <- meanheight - (t_critical_height * SEMheight)
lower_height
upper_height <- meanheight + (t_critical_height * SEMheight)
upper_height


qqPlot(height, distribution = "norm", envelope = 0.95, main = "Q-Q Plot of Tree Height", xlab = "Theoretical Normal Quantiles", ylab = "Observed Height")
skewness(height)
kurtosis(height)
shapiro.test(height)

log_height <- log(height) #log transformation of height data
qqPlot(log_height, distribution = "norm", envelope = 0.95, main = "Q-Q Plot of Log Height")
skewness(log_height)
kurtosis(log_height)
shapiro.test(log_height)

sqrt_height <- sqrt(height) #square root transformation of height data
qqPlot(sqrt_height, distribution = "norm", envelope = 0.95, main = "Q-Q Plot of Square Root Height")
skewness(sqrt_height)
kurtosis(sqrt_height)
shapiro.test(sqrt_height)


```

---

## PSet4Q2

```r
#H_0: Polarized Wings of female Heliconius butterflies had no efect on mate attraction 
#H_A: Polarized Wings of female Heliconius butterflies had an effect on mate attraction 
#Test of choice: T-test 
library(moments)
library(car)
library(ggplot2)
library(tidyverse)
library(ggpubr)

read.csv("Data/Butterfly wings.csv")
data <- read.csv("Data/Butterfly wings.csv")

#subtract Difference Depolarized - Polarized 
#Add collumn for Differences
data$Differences <- (data$Depolarized - data$Polarized)
data

hist(data$Differences, main = "Histogram of Differences in Mate Attraction", xlab = "Differences (Depolarized - Polarized)", col = "lightblue", border = "black")
qqPlot(data$Differences, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Q-Q plot of Differences in Mate Attraction", col.lines = "red")
skewness(data$Differences)  #if value for Skewness is greater than 1 or less than -1, the data is highly skewed. If the value is between -1 and -0.5 or between 0.5 and 1, the data is moderately skewed. If the value is between -0.5 and 0.5, the data is approximately symmetric.
kurtosis(data$Differences) #if value for Kurtosis is greater than 3, the data is heavy-tailed. If the value is less than 3, the data is light-tailed. If the value is approximately 3, the data is approximately normal.if 
shapiro.test(data$Differences) #If the p-value is less than 0.05, we reject the null hypothesis and conclude that the data is not normally distributed.if

#transform into long data
#take the data then gather it
datalong <- data %>%
  gather(key = "Light", value = "Number of Visits", Polarized, Depolarized)
datalong
#useggPubr for plotting
ggpaired(datalong, x = "Light", y = "Number of Visits", 
         order = c("Polarized", "Depolarized"),
         ylab  = "Number of Visits", xlab = "Light treatment")

boxplot(data$Polarized, data$Depolarized, names = c("Polarized", "Depolarized"),
        col = "skyblue", xlab = "Light treatment", ylab = "Number of Visits")

t.test(data$Polarized, data$Depolarized, paired = TRUE)
#Based on the t-test results, we can reject the null hypothesis and conclude that there is a significant difference in mate attraction between polarized and depolarized wings.ox

```

---

## PSet5Q2

```r
#ANOVA
#H_0: these conditions had no effect on plant size
#H_1: at least one condition had an effect on plant size

#load libraries
library(tidyverse)
library(car) #for QQ-plot
library(moments) #for skewness and kurtosis)
library(agricolae) #for Tukey's HSD
library(ggplot2)
library(lme4) #for random factors
library(vegan) #for perMANOVA

#Model I ANOVA
#test 

#load data
erodiumdata <- read.csv("Data/Erodium data.csv")
erodiumdata
inspect_df(erodiumdata) #check for missing data

#use data to calculate plant size as the product of height and width 
erodiumdata$productHW <- erodiumdata$Width*erodiumdata$Height
erodiumdata
inspect_df(erodiumdata)

#fit data to a linear model with Treatment = x and productHW = y
erodiumdata$Treatment <- as.factor(erodiumdata$Treatment)
lm <- lm(productHW~Treatment, data=erodiumdata) #takes the data as lm(y~x)
lm
anova(lm) #provides the ANOVA table
summary(lm) #provides extra details about the model

#check the assumptions
#recall, we want to look at the residuals to assess normality
lmres <- resid(lm)
qqPlot(lmres, "norm")
hist(lmres, xlab = "Product Height residuals", main = "Histogram")
boxplot(lmres, ylab = "Product Height residuals")
shapiro.test(lmres) #P > 0.05
skewness(lmres) 
kurtosis(lmres)
#next test equality of variances
bartlett.test(productHW ~ Treatment, data = erodiumdata) #P > 0.05
#hist(resid(lm))


#log transform the data to see if it improves normality
erodiumdata$logProduct <- log(erodiumdata$productHW) #log transform the data to see if it improves normality)
erodiumdata

#fit data to a linear model with Treatment = x and productHW = y
erodiumdata$Treatment <- as.factor(erodiumdata$Treatment)

loglm <- lm(logProduct~Treatment, data=erodiumdata) #takes the data as lm(y~x)
loglm

#check the assumptions
#recall, we want to look at the residuals to assess normality
loglmres <- resid(loglm)

qqPlot(loglmres, "norm")
hist(loglmres, xlab = "Log Product Height residuals", main = "")
boxplot(loglmres, ylab = "Log Product Height residuals")
shapiro.test(loglmres) #P > 0.05
skewness(loglmres) 
kurtosis(loglmres)
#next test equality of variances
bartlett.test(logProduct ~ Treatment, data = erodiumdata) #P > 0.05

anova(loglm) #provides the ANOVA table
summary(loglm) #provides extra details about the model


#prepare data for plotting
graphdata <- erodiumdata %>%
  group_by(Treatment) %>%
  summarize(mean=mean(productHW, na.rm=TRUE), se=sd(productHW, na.rm=TRUE)/sqrt(length(na.omit(productHW))))
graphdata
#rearrange order of the levels on the x-axis
graphdata$Treatment <- factor(graphdata$Treatment, levels = c("control", "shade", "wind", "3XH2o"))

#make a barplot with ggplot()
ggplot(data=graphdata, aes(x=Treatment, y=mean, fill=Treatment)) + 
  theme_bw() + 
  theme(axis.text.x=element_text(color="black", size=14), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=16, face="bold"), 
        axis.title.y = element_text(color="black", size=16, face="bold"),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_bar(colour="black", fill="grey", width=0.8, stat="identity") + 
  guides(fill=FALSE) + 
  ylab("Product Height x Width") +
  xlab("Treatment Type") +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), stat="identity", 
                position=position_dodge(width=0.9), width=0.1) +
  geom_text(aes(label=c("a", "bc", "b", "c"), y = mean + 1))
#note that the labels are assigned in the order of levels in "graphdata"

```

---

## PSet5Q3

```r
#H_0: the variance among the gall for mean thorax size is the same for all groups (sigma)^2_gall = 0
#H_1: the variance among the gall for mean thorax size is not the same for all groups (sigma)^2_gall > 0

aphids <- read.csv("data/Aphid data.csv")

head(aphids)
str(aphids)

#make Gall categorical
aphids$Gall <- as.factor(aphids$Gall)

model <- aov(Length ~ Gall, data = aphids)

anova_table <- summary(model)[[1]]
anova_table

#because p < 0.05, we reject the null hypothesis and conclude that there is a significant difference in mean thorax size among the different gall groups.

MS_gall <- anova_table["Gall", "Mean Sq"]   #1.7029100529
MS_within <- anova_table["Residuals", "Mean Sq"]   #0.0838095238
n <- 4 #aphids sampled per gall

variance_gall <- (MS_gall - MS_within) / n
variance_within <- MS_within
total_variance <- variance_gall + variance_within

variance_gall
variance_within
variance_gall / total_variance * 100 #genetic percent for gall
variance_within / total_variance * 100 #genetic percent for environment

#so (sigma)^2_genetic = 0.4048
#(sigma)^2_environmental = 0.0838

#(sigma)^2_total = (sigma)^2_genetic + (sigma)^2_environmental



genetic_percent <-
  variance_gall / total_variance * 100

environment_percent <-
  variance_environment / total_variance * 100

total_variance
genetic_percent
environment_percent

#conclusion: The estimated genetic variance among galls was 0.4048, while the estimated environmental variance within galls was 0.0838. 
#Thus, approximately 82.8% of the phenotypic variance in aphid thorax length was attributed to differences among genetically distinct aphid colonies, while approximately 17.2% was attributed to environmental or individual-level variation within galls.
#
```

---

## PSet5Q4 fancy

```r
#Repeated measures ANOVA
#H_0: The algal symbionts of these Polyps can not acclimate to elevated temperatures
#H_A: At least one can acclimate to elevated temperatures
#day is a repeated measure, temperature is a seperate 3 group factor, and coral is a random factor. The repeated-measures ANOVA will test for the effects of temperature, day, and their interaction on FvFm while accounting for the repeated measurements on the same coral.
#10 corals per temperature treatment, measured over 6 days. The data is balanced and there are no missing values or duplicates for the same coral on the same day.

#load libraries
library(tidyverse)
library(rstatix) #for identifying outliers
library(afex) #for ANOVA
library(ggpubr) #for plotting
library(emmeans)
library(ggplot2)
load_biostats_packages(install_missing = TRUE)

corals <- read.csv("data/Coral data.csv") %>%
  mutate(
    temperature = factor(
    temperature, 
    levels = c("ambient", "hot", "warm")
    ),
    coral = factor(coral),
    day = factor(day, levels = 1:6)
  )

# Check that the repeated-measures structure is intact
stopifnot(!anyNA(corals))

if(anyDuplicated(corals[c("coral", "day")])) {
  stop("Duplicate measurements found for the same coral on the same day.")
}

temperature_per_coral <- tapply(
  as.character(corals$temperature), corals$coral,
  function(x) length(unique(x))
  ) 

if(any(temperature_per_coral != 1)) {
  stop("At least one coral has measurements at multiple temperatures, which violates the repeated-measures assumption.")
}

cell_counts <- table(corals$temperature, corals$day)
print(cell_counts)

if(length(unique(as.vector(cell_counts))) !=1) {
  stop("The temperature-by-day groups are not balanced")
}

#Repeated-measures ANOVA:
#temperature = between-coral factor
#day = repeated within-coral factor
rm_anova <- aov_ez(
  id = "coral",
  dv = "FvFm",
  data = corals,
  between = "temperature",
  within = "day",
  anova_table = list(correction = "GG", es = "pes")
)

#Show ANOVA with Greenhouse-Geisser correction
nice(rm_anova, correction = "GG", es = "pes")

#Follow up interaction
#Compare temperatures seperately on each day
emm_temperature_by_day  <- emmeans(rm_anova, ~ temperature | day)
pairs(emm_temperature_by_day, adjust = "tukey")

#Repeated measure ANOVA with "Greenhouse-Geisser correction; needed to adjust day and temperature by day tests for potential violations of "Sphericity"
#the following comparison help show which days or temperature groups differ; the overall interaction alone does not identify every specific difference.
#the graph plots the mean FvFm for each temperature group across the 6 days, with error bars representing the standard error of the mean. The different line types and colors help distinguish between the temperature treatments.
graphdata <- corals %>%
  group_by(temperature, day) %>%
  summarise(
    mean = mean(FvFm),
    se = sd(FvFm) / sqrt(n())
  ) %>%
  mutate(day_num = as.numeric(as.character(day)))

fvfm_plot <- ggplot(
  graphdata,
  aes(
    x = day_num,
    y = mean,
    color = temperature,
    linetype = temperature, 
    group = temperature
  )
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.12,
    linewidth = 0.6
  ) +
  scale_color_manual(
    values = c(ambient = "blue", "hot" = "red", "warm" = "orange")
  ) +
  scale_linetype_manual(
    values = c(ambient = "solid", hot = "dotted", warm = "dashed")
  ) +
  scale_x_continuous(breaks = 1:6) +
  coord_cartesian(ylim = c(0.30, 0.46)) +
  labs(
    x = "Day",
    y = "PS efficiency (Fv/Fm)",
    color = "Temperature",
    linetype = "Temperature"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
fvfm_plot
#y-axis is limited to 0.30-0.46 to better visualize the differences between temperature treatments, as the FvFm values are relatively close together and this range allows for a clearer comparison of the trends over time.
#the error bars represent the standard error of the mean for each temperature group on each day, providing a visual indication of the variability and reliability of the mean estimates.

#Results
#temperature by day interaction: F(5.37, 72.49) â‰ˆ 50.38, p < 0.001 (Greenhouse-Geisser corrected). Fv/Fm patern differs over time among btemperature groups
#Day: F(2.68, 72.49) â‰ˆ 26.98, p < 0.001 Fv/Fm changes across days when averaged over temperatures
#temperature averaged across days: F(2, 27) â‰ˆ 14.79, p < 0.001, this overall comaprison is not significant at a 0.05 level
#Interpretation: the key result is the Interaction, where the temperatures do not show the same pattern throughout the 6 days. This supports a temperature dependent change over time. 
#Conclusion: Repeated Measure ANOVA found significant affects of temperature, day, and their interaction (all p < 0.05). interaction means Fv/Fm changed differently over time across temperatures: initial decreases followed by increases over time, suggestingng they are acclimating! 
```

---

## PSet5Q5

```r
diversity_data <- read.csv("data/Diversity data.csv")
diversity_data

#location is treated as category
diversity_data$Location <- as.factor(diversity_data$Location)

#diversity remains numeric
str(diversity_data)

#sample variance within each Location
variance_by_location <- tapply(
  diversity_data$Diversity,
  diversity_data$Location,
  var
)
variance_by_location

inspect_df(diversity_data)

#diversity_data$Diversity <- as.factor(diversity_data$Diversity)
lm_model <- lm(Diversity ~ Location, data = diversity_data)
lm_model

anova(lm_model)
summary(lm_model)

lm_model_res <- resid(lm_model)

qqPlot(lm_model_res, "norm")
hist(lm_model_res, xlab = "Diversity residuals", main = "")
boxplot(lm_model_res, ylab = "Diversity residuals")
shapiro.test(lm_model_res) #P > 0.05
skewness(lm_model_res) 
kurtosis(lm_model_res)
#next test equality of variances
bartlett.test(lm_model_res ~ Location, data = diversity_data) #P > 0.05

diversity_data$log_diversity_data <- log(diversity_data$Diversity)

log_lm_model <- lm(log_diversity_data ~ Location, data = diversity_data)
log_lm_model_res <- resid(log_lm_model)
summary(log_lm_model)

qqPlot(log_lm_model_res, "norm")
hist(log_lm_model_res, xlab = "Log Diversity residuals", main = "")
boxplot(log_lm_model_res, ylab = "Log Diversity residuals")
shapiro.test(log_lm_model_res) #P > 0.05
skewness(log_lm_model_res) 
kurtosis(log_lm_model_res)
#next test equality of variances
bartlett.test(log_lm_model_res ~ Location, data = diversity_data) #P > 0.05

diversity_data$sqrt_diversity_data <- sqrt(diversity_data$Diversity)
sqrt_lm_model <- lm(sqrt_diversity_data ~ Location, data = diversity_data)
sqrt_lm_model_res <- resid(sqrt_lm_model)
summary(sqrt_lm_model)

qqPlot(sqrt_lm_model_res, "norm")
hist(sqrt_lm_model_res, xlab = "Sqrt Diversity residuals", main = "")
boxplot(sqrt_lm_model_res, ylab = "Sqrt Diversity residuals")
shapiro.test(sqrt_lm_model_res) #P > 0.05
skewness(sqrt_lm_model_res) 
kurtosis(sqrt_lm_model_res)
#next test equality of variances
bartlett.test(sqrt_lm_model_res ~ Location, data = diversity_data) #P > 0.05

kruskal.test(Diversity ~ Location, data = diversity_data)
bartlett.test(Diversity ~ Location, data = diversity_data)

```

---

## PSet6Q1

```r
load_biostats_packages()
#load libraries
library(tidyverse)
library(car)
library(ggplot2)
library(GGally)
library(psych)
library(moments)
#H_0: Beak width and beak depth have no linear association in the population; the correlation is zero (Ï = 0)
#H_A: Beak width and beak depth have a linear association; the correlation is not zero (Ï â‰  0)
penguins <- read.csv("data/penguin data - Sheet1.csv")
View(penguins)

#penguins$Beak.width <- as.numeric(penguins$Beak.width)
qqPlot(penguins$Beak.width, "norm",ylab = "", main = "Q-Q plot of Beak width")
skewness(penguins$Beak.width) #near 0 = more symmetrical, + = longer right-tail, - = longer left-tail 
kurtosis(penguins$Beak.width) #normal distribution ~ 3, >3 = more peaked, <3 = flatter
shapiro.test(penguins$Beak.width) #p-value > 0.05 = no strong evidence against normality, p-value <= 0.05 = evidence against normality
hist(penguins$Beak.width, main = "Histogram of Beak width", xlab = "Beak width (mm)", ylab = "Frequency", col = "lightblue", border = "black")

#penguins$Beak.width <- as.numeric(penguins$Beak.width)
qqPlot(penguins$Beak.depth, "norm", ylab = "", main = "Q-Q plot of Beak depth")
skewness(penguins$Beak.depth) #near 0 = more symmetrical, + = longer right-tail, - = longer left-tail 
kurtosis(penguins$Beak.depth) #normal distribution ~ 3, >3 = more peaked, <3 = flatter
shapiro.test(penguins$Beak.depth) #p-value > 0.05 = no strong evidence against normality, p-value <= 0.05 = evidence against normality
hist(penguins$Beak.depth, main = "Histogram of Beak depth", xlab = "Beak depth (mm)", ylab = "Frequency", col = "lightblue", border = "black")

#scatterplot the relationship between beak width and beak depth
plot(penguins$Beak.width ~ penguins$Beak.depth, data = penguins,
     xlab = "Beak depth (mm)",
     ylab = "Beak width (mm)",
     main = "Beak width by Beak depth in African penguins",
     pch = 19)
abline(lm(penguins$Beak.width ~ penguins$Beak.depth, data = penguins), col = "red")

cor.test(penguins$Beak.depth, penguins$Beak.width, method = "pearson")

#penguins_long <- penguins %>%
#  gather(key = "Trait", value = "Size", Beak.width, Beak.depth)
#penguins_long 
#view(penguins_long)
#bartlett.test(Size ~ Trait, data = penguins_long)


#plot(penguins$Beak.width, penguins$Beak.depth,
#     xlab = "Beak width (mm)",
#     ylab = "Beak depth (mm)",
#     main = "Beak width and depth in African penguins",
#     pch = 19)
#abline(lm(penguins$Beak.depth ~ penguins$Beak.width), col = "blue", lwd = 2)

#model <- lm(penguins$Beak.depth ~ penguins$Beak.width)

#plot(model, which = 1)  
#Plot 1: Residuals vs. Fitted. residual is (depth y_i) - (depth y-hat), looking for residuals scattered random around 0, curve ~suggest relationship isnâ€™t linear; funnel shape ~suggest residual spread changes across fitted values
#plot(model, which = 2)  
#Plot 2: Normal Q-Q. compares residuals with what weâ€™d expect if they were normally distributed. Points roughly following line suggest normality is reasonable, visually complements Shapiroâ€“Wilk test ran

#shapiro.test(resid(model))


#penguins$differences <- penguins$Beak.depth - penguins$Beak.width
#var(penguins$differences)
#qqPlot(penguins$differences)
#shapiro.test(penguins$differences)

```

---

## PSet6Q2

```r
#clear the environment
rm(list=ls())

#load libraries
library(tidyverse)
library(car)
library(ggplot2)
library(GGally)
library(psych)

edata <- read.csv("data/Environmental data.csv")
view(edata)
str(edata)

data_subset <-  edata[, c("PctCover", "Elevation", "TPI", "TWI", "MaxTemp")]

#create the scatter plot matrix
pairs(data_subset, 
      main = "SPLOM",
      pch = 20,                                     
      cex = 1.5,                                    
      col = "#69b3a2"                               
)
ggpairs(data_subset, 
        title = "")

#yet another method using pair.panels() from library(psych)
#this gives histograms and ellipses
pairs.panels(data_subset,
             method = "pearson", # Correlation method (e.g., pearson, spearman)
             hist.col = "#00AFBB", # Color for histograms
             density = TRUE,      # Add density plots
             ellipses = TRUE      # Add correlation ellipses
)

#create the individual Pearson's correlation tests
#notice that the results match those shown above
cor.test(edata$PctCover, edata$Elevation, method="pearson", na.rm=TRUE)
cor.test(edata$PctCover, edata$TPI, method="pearson", na.rm=TRUE)
cor.test(edata$PctCover, edata$TWI, method="pearson", na.rm=TRUE)
cor.test(edata$PctCover, edata$MaxTemp, method="pearson", na.rm=TRUE)
cor.test(edata$Elevation, edata$TPI, method="pearson", na.rm=TRUE)
cor.test(edata$Elevation, edata$TWI, method="pearson", na.rm=TRUE)
cor.test(edata$Elevation, edata$MaxTemp, method="pearson", na.rm=TRUE) #***
cor.test(edata$TPI, edata$TWI, method="pearson", na.rm=TRUE)
cor.test(edata$TPI, edata$MaxTemp, method="pearson", na.rm=TRUE)
cor.test(edata$TWI, edata$MaxTemp, method="pearson", na.rm=TRUE)


```

---

## PSet6Q4 3 predictors seperately

```r
load_biostats_packages()
data <- read.csv("data/Krat data.csv")
View(data)
inspect_df(data)

#fit three separate models
#model for Shrub Cover
shrub_model <- lm(RatDnsty ~ ShrbCvr, data = data,
                  na.action = "na.omit")

#model for Seed Production
seed_model <- lm(RatDnsty ~ SeedPrd, data = data,
                 na.action = "na.omit")

#model for Snake Density
snake_model <- lm(RatDnsty ~ SnkDnsty, data = data,
                  na.action = "na.omit")

#check Shrub Cover model for assumptions 
shrub_modelres <- resid(shrub_model)
qqp(shrub_modelres, "norm", main = "qqPlot Shrub Cover Model Residuals")
hist(shrub_modelres, main = "Histogram Shrub Cover Model Residuals", xlab = "Residuals")
shapiro.test(shrub_modelres)

#plot(shrub_modelres ~ fitted(shrub_model),
#     ylab = "Residuals",
#     xlab = "Predicted rat density")
#abline(h = 0, col = "red", lty = 2)

#make a scatterplot with trendline
ggplot(data, aes(x=RatDnsty, y=ShrbCvr)) +
  geom_point(size = 2) +
  xlab("Rat Density") +
  ylab("Shrub Cover") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_smooth(method = "lm") #adds regression line with 95% confidence intervals

#check Seed Production model for assumptions 
seed_modelres <- resid(seed_model)
qqp(seed_modelres, "norm", main = "qqPlot Seed Production Model Residuals")
hist(seed_modelres, main = "Histogram Seed Production Model Residuals", xlab = "Residuals")
shapiro.test(seed_modelres)

#plot(seed_modelres ~ fitted(seed_model),
#     ylab = "Residuals",
#     xlab = "Predicted rat density")
#abline(h = 0, col = "red", lty = 2)

ggplot(data, aes(x=RatDnsty, y=SeedPrd)) +
  geom_point(size = 2) +
  xlab("Rat Density") +
  ylab("Seed Production") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_smooth(method = "lm") #adds regression line with 95% confidence intervals

#check Snake Density model for assumptions 
snake_modelres <- resid(snake_model)
qqp(snake_modelres, "norm", main = "qqPlot Snake Density Model Residuals")
hist(snake_modelres, main = "Histogram Snake Density Model Residuals", xlab = "Residuals")
shapiro.test(snake_modelres)

#plot(snake_modelres ~ fitted(snake_model),
#     ylab = "Residuals",
#     xlab = "Predicted rat density")
#abline(h = 0, col = "red", lty = 2)

ggplot(data, aes(x=RatDnsty, y=SnkDnsty)) +
  geom_point(size = 2) +
  xlab("Rat Density") +
  ylab("Snake Density") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_smooth(method = "lm") #adds regression line with 95% confidence intervals

#assumptions for models relatively satisfied, examine results
summary(shrub_model)
anova(shrub_model)
summary(seed_model)
anova(seed_model)
summary(snake_model)
anova(snake_model)

assumption_report(data = shrub_modelres)
assumption_report(data = seed_modelres)
assumption_report(data = snake_modelres)

```

---

## PSet6Q4 3 predictors together Regression

```r
#This asks whether each factor is related to rat density while accounting for the other two factors

load_biostats_packages()
data <- read.csv("data/Krat data.csv")
view(data)
analysis_guide(data, outcome = "RatDnsty", predictors = c("ShrbCvr", "SeedPrd", "SnkDnsty"))
#Fit one model with all three predictors
datalin <- lm(
  RatDnsty ~ ShrbCvr + SeedPrd + SnkDnsty, #~ means â€œas a function ofâ€
  data = data,
  na.action = "na.omit"
)

#Save this model's residuals
datalinres <- resid(datalin) #resid(datalin) extracts the differences between the observed rat densities and the values predicted by this model

#Check whether the residuals are reasonably normal
qqp(datalinres, "norm")
hist(datalinres)
shapiro.test(datalinres)

#Check residual spread against predicted rat density
#plot(datalinres ~ fitted(datalin)) plots residuals vertically against predicted values horizontally. Look for a roughly even cloud around zero, without a curve or funnel shape
plot(datalinres ~ fitted(datalin), #fitted(datalin) extracts the modelâ€™s predicted rat densities
     ylab = "Residuals",
     xlab = "Predicted rat density")
abline(h = 0, col = "red", lty = 2)

#Examine the model's results
summary(datalin)

termplot(datalin, se = TRUE, partial.resid = FALSE)

```

---

## PSet6Q4

```r
load_biostats_packages()
data <- read.csv("data/Krat data.csv")
view(data)

#first fit a linear model
datalin <- lm(RatDnsty ~ ShrbCvr + SeedPrd + SnkDnsty, data=data, na.action="na.omit")
mfdatalin <- model.frame(datalin)

datalinres <- resid(datalin)
#termplot(datalin, se = TRUE, partial.resid = FALSE)
#plot model
plot(RatDnsty ~ ShrbCvr, data = data,
     xlab = "Shrub cover", ylab = "Rat density")
abline(lm(data$RatDnsty ~ data$ShrbCvr, data = data), col = "red")
plot(RatDnsty ~ SeedPrd, data = data,
     xlab = "Seed production", ylab = "Rat density")
abline(lm(data$RatDnsty ~ data$SeedPrd, data = data), col = "red")
plot(RatDnsty ~ SnkDnsty, data = data,
     xlab = "Snake density", ylab = "Rat density")
abline(lm(data$RatDnsty ~ data$SnkDnsty, data = data), col = "red")

#check the assumptions
#recall, we want to look at the residuals
#first test normality
datalinres <- residuals(datalin)
#qqPlot(datalinres, main = "Q-Q plot of residuals", ylab = "Residuals")

#check the assumptions
#recall, we want to look at the residuals
#first test normality
#model1res <- resid(model1)
#qqp(model1res, "norm")
#hist(model1res)
#shapiro.test(model1res) #P < 0.05. It's a bit non-normal
#test for homogeneity of variance
#call the fitted y values
#fitted(model1) #this is the predicted value of y, for each value of x
#plot the fitted (predicted) values against the residuals
#plot(model1res~fitted(model1), ylab = "Residuals", xlab = "Predicted y") #it looks a bit cone shaped


shrub_model <- lm(RatDnsty ~ ShrbCvr, data = data)
qqp(residuals(shrub_model), "norm")
hist(residuals(shrub_model))
shapiro.test(residuals(shrub_model)) #P > 0.05. This is better
fitted(shrub_model)
plot(residuals(shrub_model) ~ fitted(shrub_model), ylab = "Residual", xlab = "Predicted y")

ggplot(data, aes(x = ShrbCvr, y = RatDnsty)) +
  geom_point() +
  xlab("Shrub cover") +
  ylab("Rat density") +
  geom_smooth(method = "lm", formula = y ~ x, 
              se = TRUE, level = 0.95) +
  labs(x = "Shrub cover", y = "Rat density")

seed_model <- lm(RatDnsty ~ SeedPrd, data = data)
qqp(residuals(seed_model), "norm")
hist(residuals(seed_model))
shapiro.test(residuals(seed_model)) #P > 0.05. This is better
fitted(seed_model)
plot(residuals(seed_model) ~ fitted(seed_model), ylab = "Residual", xlab = "Predicted y")

ggplot(data, aes(x = SeedPrd, y = RatDnsty)) +
  geom_point() +
  xlab("Seed production") +
  ylab("Rat density") +
  geom_smooth(method = "lm", formula = y ~ x, 
              se = TRUE, level = 0.95) +
  labs(x = "Seed production", y = "Rat density")

snake_model <- lm(RatDnsty ~ SnkDnsty, data = data)
qqp(residuals(snake_model), "norm")
hist(residuals(snake_model))
shapiro.test(residuals(snake_model)) #P > 0.05. This is better
fitted(snake_model)
plot(residuals(snake_model) ~ fitted(snake_model), ylab = "Residual", xlab = "Predicted y")

ggplot(data, aes(x = SnkDnsty, y = RatDnsty)) +
  geom_point() +
  xlab("Snake density") +
  ylab("Rat density") +
  geom_smooth(method = "lm", formula = y ~ x, 
              se = TRUE, level = 0.95) +
  labs(x = "Snake density", y = "Rat density")


#plot(datalinres ~ data$ShrbCvr, xlab = "Shrub cover", ylab = "Residuals")


#try a log transformation
data$logRatDnsty <- log(data$RatDnsty)
data$logShrbCvr <- log(data$ShrbCvr)
data$logSeedPrd <- log(data$SeedPrd)
data$logSnkDnsty <- log(data$SnkDnsty)
data

#create a new linear model with the log data
logmodel<-lm(logNumberSpp~logBiomass, data=riffledata, na.action="na.omit")

#check the assumptions of the log model
#test normality
logmodelres <- resid(logmodel)
qqp(logmodelres, "norm")
hist(logmodelres)
shapiro.test(logmodelres) #P > 0.05. This is better
#test for homogeneity of variance
# call the fitted y values
fitted(logmodel) 
#plot the fitted values against the residuals
plot(logmodelres~fitted(logmodel), ylab = "Residuals", xlab = "Predicted log y") #it looks better

#let's perfrom regression with the transformed data
summary(logmodel) #notice it also does a t-test (but t is calculated differently than in correlation)
#make an ANOVA table
anova(logmodel)


#let's perfrom regression with the transformed data
summary(logmodel) #notice it also does a t-test (but t is calculated differently than in correlation)
#make an ANOVA table
anova(logmodel)

#make a scatterplot with trendline
ggplot(riffledata, aes(x=logBiomass, y=logNumberSpp)) +
  geom_point(size = 2) +
  xlab("Log (Biomass)") +
  ylab("Log (Species Richness)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_smooth(method = "lm") #adds regression line with 95% confidence intervals

```

---

## PSet6Q5

```r
load_biostats_packages()
data <- read.csv("data/Hawk data.csv")
view(data)
str(data)
inspect_df(data)

#check assumptions of normality
hist(data$Jawforce)
shapiro.test(data$Jawforce)
qqPlot(data$Jawforce)
assumption_report(data = data)

hist(data$Footforce)
shapiro.test(data$Footforce)
qqPlot(data$Footforce)

cor.test(data$Jawforce, data$Footforce, method = "pearson")

#make a scatterplot
ggplot(data, aes(x=Jawforce, y=Footforce)) +
  geom_point(size = 2) +
  xlab("Footforce") +
  ylab("Jawforce") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
abline(lm(Footforce ~ Jawforce, data = data), col = "red")

#scatterplot the relationship between jawforce and footforce 
plot(Jawforce ~ Footforce, data = data,
     xlab = "Footforce",
     ylab = "Jawforce",
     main = "Footforce by Jawforce Hawks",
     pch = 19)
abline(lm(Jawforce ~ Footforce, data = data), col = "red")

correlation_workflow(data, x = "Jawforce", y = "Footforce", method = "pearson", plots = TRUE)

```

---

## QSet5Q4

```r
#Repeated measures ANOVA
#H_0: The algal symbionts of these Polyps can not acclimate to elevated temperatures
#H_A: At least one can acclimate to elevated temperatures

#load libraries
library(tidyverse)
library(rstatix) #for identifying outliers
library(afex) #for ANOVA
library(ggpubr) #for plotting
library(emmeans)
library(ggplot2)
load_biostats_packages(install_missing = TRUE)

corals <- read.csv("data/Coral data.csv")
inspect_df(corals)
corals

#check model assumptions
qqPlot(corals$FvFm, ylab = "FvFm")
shapiro.test(corals$FvFm) #P = 0.1, data has normal distribution
bartlett.test(FvFm ~ day, data = corals) #P = 0.1, data has equal variance
corals %>%
  group_by(day) %>%
  identify_outliers(FvFm) #no outliers

#perform repeated measures ANOVA
#build linear model
rm_anova <- aov_ez(data = corals, dv = "FvFm", id = "coral", within = "day")
rm_res <- resid(rm_anova, append = TRUE)
hist(rm_res$.residuals, breaks = 20, xlab = "Residuals FvFm", ylab = "Frequency", main = "")
qqPlot(rm_res$.residuals, ylab = "FvFm residuals")
shapiro.test(rm_res$.residuals) #P = 0.1, data has normal distribution
skewness(rm_res$.residuals)
kurtosis(rm_res$.residuals)
bartlett.test(FvFm ~ day, data = corals)
bartlett.test(rm_res$.residuals ~ corals$day, data = corals)
boxplot(rm_res$.residuals, ylab = "FvFm residuals")
var(rm_res$.residuals)
var(corals$FvFm)
inspect_df(rm_res)

#view the ANOVA & HSD results
summary(rm_anova, append = TRUE)
em_results <- emmeans(rm_anova, ~ day)
pairs(em_results, adjust = "tukey") #Tukey's HSD

#make a boxplot
bxp <- ggboxplot(corals, x = "day", y = "FvFm", add = "point", notch = TRUE)
bxp

#try a bar plot
#calculate means and variance
graphdata <- corals %>%
  group_by(day) %>%
  summarize(mean = mean(FvFm, na.rm = TRUE), se = sd(FvFm, na.rm = TRUE) / sqrt(length(na.omit(FvFm))))
graphdata

#ggplot
ggplot(data = graphdata, aes(x = Treatment, y = mean, fill = Treatment)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 16), 
        axis.text.y = element_text(face = "bold", color = "black", size = 13), 
        axis.title.x = element_text(color = "black", size = 14, face = "bold"), 
        axis.title.y = element_text(color = "black", size = 14, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_bar(colour = "black", fill = "DodgerBlue", width = 0.7, stat = "identity") + 
  guides(fill = FALSE) + 
  ylab("Average FvFm") +
  xlab("") +
  scale_x_discrete(labels = c("baseline" = "day 1", "1" = "day 2", "2" = "day 3", "3" = "day 4", "4" = "day 5", "5" = "day 6")) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), stat = "identity", 
                position = position_dodge(width = 0.9), width = 0.1) +
  coord_cartesian(ylim = c(75, 100))


```

---

## R script template

```r
#R script template

#Getting a new project going
#It's imperative to have an organized directory to make R work smoothly

#FIRST
#Create three new folders within your project's directory
#This can be done outside of R studio, or in the Files window at lower right
#Name the folders:
#Data
#Script
#Output

#You can copy the first two scripts into the Script folder outside of R studio

#The first line of an R Script should be some kind of title
#Followed by the author's name and date
#For example:

#New R Script
#Created by you
#Date

#Next provide a brief description of what the code is going to do
#We'll work with the Mytilus byssal thread break force data here for practice

#Often times, a script will begin by sweeping the environment clean
rm(list=ls())
#This can also be done by clicking the broom above the environment window

#SECOND
#Load libraries
#Libraries are software packages
#Click on the Packages tab at lower right to see a list of them
#You can install them by clicking on the box to the left
#Or by using the function install.packages()
#For example:
install.packages('tidyverse')
#Notice that the package is being installed in the console
#You can usually ignore this until it's done 
#Now load the package with the library() function
library(tidyverse)
#The tidyverse package has many useful functions for data cleaning
#We'll use lots of other packages eventually
#For now the tidyverse will help us bring an Excel spreadsheet into R

#THIRD
#Load data
#Let's start with loading the spreadsheet created in Excel
#We are going to make it a new object
data <- read.csv("Data/Byssal threads(sheet1).csv")
#read.csv calls the spreadsheet
#The directory must be specified exactly or it won't work
#You can find your working directory with getwd()
getwd()
#View the data in a new tab with View()
View(data)
#Or just run it in the console
data
#Look at the first few rows only
head(data)
#Look at the data's structure
str(data)

#FOURTH
#We can calculate summary statistics by site
#We need the %>% from the tidyverse (this is called a pipe)
#Pristine site
prmeanstr <- data %>%
  filter(Site=="Pristine") %>%
  summarise(mean(Strength)) #returns the mean strength at the pristine site
prstdstr <- data %>%
  filter(Site=="Pristine") %>%
  summarise(sd(Strength)) #returns the stdev strength at the pristine site
#Polluted site
plmeanstr <- data %>%
  filter(Site=="Polluted") %>%
  summarise(mean(Strength)) #returns the mean strength at the polluted site
plstdstr <- data %>%
  filter(Site=="Polluted") %>%
  summarise(sd(Strength)) #returns the stdev strength at the polluted site

#We can make a simple graph
#first set means as numeric
prmeanstr <- as.numeric(prmeanstr)
plmeanstr <- as.numeric(plmeanstr)
#concatenate them into a variable y
y <- c(prmeanstr, plmeanstr)
#name the levels of x
names <- matrix(c("Pristine", "Polluted"))
#make a data frame
graphdata <- data.frame(Site = names, Strength = y)
#use barplot() for graphing
barplot(height = graphdata$Strength, names.arg = graphdata$Site, 
        ylim = c(0,1000), xlab = "Site", ylab = "Mean Strength (g)")

#Add error bars with the standard deviations using the arrows() function
#first make stdev numeric
prstdstr <- as.numeric(prstdstr)
plstdstr <- as.numeric(plstdstr)
#concatenate them into an object
stdev <- c(prstdstr, plstdstr)
#redefine the plot with greater y axis limit that fits bars
barploter <- barplot(height = graphdata$Strength, names.arg = graphdata$Site, 
                     ylim = c(0,1500), xlab = "Site", ylab = "Mean Strength (g)")
#add bars
arrows(x0 = barploter, y0 = y - stdev,
       x1 = barploter, y1 = y + stdev,
       angle = 90, code = 3, length = 0.05)

#calculate confidence limits
#first calculate summary stats for the entire data set
ybar <- mean(data$Strength)
SEMy <- sd(data$Strength)/sqrt(count(data))
#find the critical value of t with n-1 df
qt(p = 0.025, df = 24, lower.tail = FALSE)
lowerCL <- ybar-(SEMy*2.063899)
lowerCL
upperCL <- ybar+(SEMy*2.063899)
upperCL

#make a histogram w/ the confidence limits
hist(data$Strength, breaks = 8, main = "", xlab= "Strength")
abline(v=lowerCL, col="blue")
abline(v=upperCL, col="blue")

```

---

## Regression

```r
#Linear regression

#clear the environment
rm(list=ls())

#load libraries
library(tidyverse)
library(car) #for QQ-plot
library(ggplot2)

#let's analyze the species richness~biomass relationship

#load data
riffledata <- read.csv("Data/Stream data.csv")
riffledata
view(riffledata)
#first fit a linear model
model1 <- lm(NumberSpp~Biomass, data=riffledata, na.action="na.omit")
#plot model
plot(NumberSpp~Biomass, data=riffledata, col="blue", ylab="NumberSpp (mg)", xlab="Biomass (cm)")
abline(model1, col="blue") #adds fit from the model

#check the assumptions
#recall, we want to look at the residuals
#first test normality
model1res <- resid(model1)
qqp(model1res, "norm")
hist(model1res)
shapiro.test(model1res) #P < 0.05. It's a bit non-normal
#test for homogeneity of variance
#call the fitted y values
fitted(model1) #this is the predicted value of y, for each value of x
#plot the fitted (predicted) values against the residuals
plot(model1res~fitted(model1), ylab = "Residuals", xlab = "Predicted y") #it looks a bit cone shaped

#try a log transformation
riffledata$logBiomass<-log(riffledata$Biomass)
riffledata$logNumberSpp<-log(riffledata$NumberSpp)
riffledata

#create a new linear model with the log data
logmodel<-lm(logNumberSpp~logBiomass, data=riffledata, na.action="na.omit")

#check the assumptions of the log model
#test normality
logmodelres <- resid(logmodel)
qqp(logmodelres, "norm")
hist(logmodelres)
shapiro.test(logmodelres) #P > 0.05. This is better
#test for homogeneity of variance
# call the fitted y values
fitted(logmodel) 
#plot the fitted values against the residuals
plot(logmodelres~fitted(logmodel), ylab = "Residuals", xlab = "Predicted log y") #it looks better

#let's perfrom regression with the transformed data
summary(logmodel) #notice it also does a t-test (but t is calculated differently than in correlation)
#make an ANOVA table
anova(logmodel)

#make a scatterplot with trendline
ggplot(riffledata, aes(x=logBiomass, y=logNumberSpp)) +
  geom_point(size = 2) +
  xlab("Log (Biomass)") +
  ylab("Log (Species Richness)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
        geom_smooth(method = "lm") #adds regression line with 95% confidence intervals

```

---

## Repeated measures ANOVA

```r
#Repeated measures ANOVA

#load libraries
library(tidyverse)
library(rstatix) #for identifying outliers
library(afex) #for ANOVA
library(ggpubr) #for plotting
library(emmeans)
library(ggplot2)

#load data
biofeedback <- read.csv("Data/Stress data.csv")
biofeedback
inspect

#check model assumptions
qqPlot(biofeedback$AverageBPM, ylab = "Average BPM")
shapiro.test(biofeedback$AverageBPM) #P = 0.1, data has normal distribution
biofeedback %>%
  group_by(Test) %>%
  identify_outliers(AverageBPM) #no outliers

#perform repeated measures ANOVA
#build linear model
rm_anova <- aov_ez(data = biofeedback, dv = "AverageBPM", id = "Subject", within = "Test")
#view the ANOVA & HSD results
summary(rm_anova)
em_results <- emmeans(rm_anova, ~ Test)
pairs(em_results, adjust = "tukey") #Tukey's HSD

#make a boxplot
bxp <- ggboxplot(biofeedback, x = "Test", y = "AverageBPM", add = "point", notch = FALSE)
bxp

#try a bar plot
#calculate means and variance
graphdata <- biofeedback %>%
  group_by(Test) %>%
  summarize(mean=mean(AverageBPM, na.rm=TRUE), se=sd(AverageBPM, na.rm=TRUE)/sqrt(length(na.omit(AverageBPM))))
graphdata

#ggplot
ggplot(data=graphdata, aes(x=Test, y=mean, fill=Test)) + 
  theme_bw() + 
  theme(axis.text.x=element_text(face="bold", color="black", size=16), 
        axis.text.y=element_text(face="bold", color="black", size=13), 
        axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_bar(colour="black", fill="DodgerBlue", width=0.7, stat="identity") + 
  guides(fill=FALSE) + 
  ylab("Average Heart Rate (BPM)") +
  xlab("") +
  #scale_x_discrete(labels=c("Baseline" = "Baseline", "1" = " Test 1", "2" = "Test 2", "3" = "Test 3")) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), stat="identity", 
                position=position_dodge(width=0.9), width=0.1) +
  coord_cartesian(ylim = c(75, 100))

```

---

## reset_workspace

```r
reset_workspace <- function() {
  rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  gc(verbose = FALSE)
  cat("\014") # Clears the console screen
  message("Global environment cleared and memory freed.")
}
```

---

## RScript Mean Thread.count

```r
#Imports Spreadsheet
data <- read.csv("Data/Byssal threads(sheet1).csv")

#Opens Spreadsheet In New Tab
View(data)

#Calculates Mean Thread.count for Pristine Site
prmeanthr <- data %>%
  filter(Site=="Pristine") %>%
  summarise(mean(Thread.count))

#Calculates Thread.count Standard Deviation for Pristine Site
prstdthr <- data %>%
  filter(Site=="Pristine") %>%
  summarise(sd(Thread.count))

#Calculates Mean Thread.count for Polluted Site 
pomeanthr <- data %>%
  filter(Site=="Polluted") %>%
  summarise(mean(Thread.count)) 

#Calculates Thread.count Standard Deviation for Polluted Site
postdthr <- data %>%
  filter(Site=="Polluted") %>%
  summarise(sd(Thread.count))
#Summary Function Above Generated Summary of Respective Data in a "mini-Table"

#Converts "mini-Table" Summary into Numeric Value
prmeanthr <- as.numeric(prmeanthr)
pomeanthr <- as.numeric(pomeanthr)

#Combines two values into single vector
y <- c(prmeanthr, pomeanthr)

#Creates site labels under bar graph - by making them into single line matrix 
names <- matrix(c("Pristine", "Polluted"))

#Creates data table that'll be used for the graph
graphdata <- data.frame(Site = names, Thread.count = y)

#Creates the actual bar graph 
barplot(height = graphdata$Thread.count, names.arg = graphdata$Site, 
        
        #sets y-axis from 0 to 60        
        ylim = c(0,60), xlab = "Site", ylab = "Mean Thread.count")

#Turns Pristine & Polluted SD from mini-table into a number
prstdthr <- as.numeric(prstdthr)
postdthr <- as.numeric(postdthr)

#Combines the two numbers into single vector
stdev <- c(prstdthr, postdthr)

#Redraws graph & saves bar position into barploter
barploter <- barplot(height = graphdata$Thread.count, names.arg = graphdata$Site, 
                     ylim = c(0,60), xlab = "Site", ylab = "Mean Thread.count")

#Draws Error bars 
arrows(x0 = barploter, y0 = y - stdev,
       x1 = barploter, y1 = y + stdev,
       angle = 90, code = 3, length = 0.05)
```

---

## snips

```r
krats <- read.csv("data/Krat data.csv")
view(krats)

vari <- c("RatDnsty", "ShrbCvr", "SeedPrd", "SnkDnsty")
dat <- krats[complete.cases(krats[vari]), vari] #complete.cases(...) checks each row, returns TRUE if it has a value for every selected column, FALSE if anyis missing
scatterplotMatrix(dat)

pairs <- combn(vari, 2, simplify = FALSE)

```

---

## tricks

```r
#Run these two together to get regexplain powers
install.packages("remotes")
remotes::install_github("gadenbuie/regexplain")

#
install.packages("easystats")
library(easystats)

install.packages("gtsummary")
library(gtsummary)

install.packages("esquisse")
# Open via RStudio Addins menu -> "ggplot2 builder"

#Excel files integration easily 
install.packages("readxl")
library(readxl)

# Descriptive baseline table split by treatment group:
tbl_summary(data, by = treatment) |> add_p()

install.packages("rex")
library(rex)
install.packages("unglue")
library(unglue)
install.packages("tidyverse")
library(tidyverse)
install.packages("car")
library(car)
install.packages("moments")
library(moments)
install.packages("pubr")
library(pubr)
install.packages("ggplot2")
library(ggplot2)

inspect_df <- function(df, digits = 2) {
  # --- 1. Dataset Overview ---
  cat("========================================================\n")
  cat("                  DATASET OVERVIEW                      \n")
  cat("========================================================\n")
  cat(sprintf("Dimensions : %d rows x %d columns\n", nrow(df), ncol(df)))
  total_cells <- prod(dim(df))
  total_na <- sum(is.na(df))
  cat(sprintf("Total NA   : %d (%.1f%% of all cells)\n", total_na, (total_na / total_cells) * 100))
  cat(sprintf("Completeness: %d complete cases (%.1f%%)\n\n", 
              sum(complete.cases(df)), 
              (sum(complete.cases(df)) / nrow(df)) * 100))
  
  # Helper for statistical mode
  calc_mode <- function(x) {
    clean_x <- x[!is.na(x)]
    if (length(clean_x) == 0) return(NA)
    ux <- unique(clean_x)
    tab <- tabulate(match(clean_x, ux))
    ux[tab == max(tab)][1] # returns first mode if multimodal
  }
  
  # --- 2. Numeric Summary Table ---
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) > 0) {
    cat("--- Numeric Variables (Parametric & Non-Parametric) ---\n")
    num_list <- lapply(num_cols, function(col) {
      vals <- df[[col]]
      clean <- vals[!is.na(vals)]
      n_clean <- length(clean)
      
      if (n_clean > 1) {
        m    <- mean(clean)
        med  <- median(clean)
        mod  <- calc_mode(clean)
        s_d  <- sd(clean)
        vari <- var(clean)
        sem  <- s_d / sqrt(n_clean)
        iqr_val <- IQR(clean)
        mini <- min(clean)
        maxi <- max(clean)
      } else {
        m <- med <- mod <- s_d <- vari <- sem <- iqr_val <- mini <- maxi <- NA
      }
      
      data.frame(
        Variable = col,
        N = n_clean,
        NA_Count = sum(is.na(vals)),
        Mean = round(m, digits),
        SD = round(s_d, digits),
        SEM = round(sem, digits),
        Var = round(vari, digits),
        Median = round(med, digits),
        IQR = round(iqr_val, digits),
        Min = round(mini, digits),
        Max = round(maxi, digits),
        Mode = round(mod, digits),
        stringsAsFactors = FALSE
      )
    })
    print(do.call(rbind, num_list), row.names = FALSE)
    cat("\n")
  }
  
  # --- 3. Categorical / Factor Summary Table ---
  cat_cols <- names(df)[!sapply(df, is.numeric)]
  if (length(cat_cols) > 0) {
    cat("--- Categorical & Factor Variables ---\n")
    cat_list <- lapply(cat_cols, function(col) {
      vals <- df[[col]]
      clean <- vals[!is.na(vals)]
      top_val <- if (length(clean) > 0) names(sort(table(clean), decreasing = TRUE))[1] else NA
      top_freq <- if (length(clean) > 0) max(table(clean)) else 0
      
      data.frame(
        Variable = col,
        Type = class(vals)[1],
        Unique = length(unique(clean)),
        NA_Count = sum(is.na(vals)),
        Top_Category = as.character(top_val),
        Top_Freq = sprintf("%d (%.1f%%)", top_freq, (top_freq / max(1, length(clean))) * 100),
        stringsAsFactors = FALSE
      )
    })
    print(do.call(rbind, cat_list), row.names = FALSE)
  }
  cat("========================================================\n")
}
```

---

## Two sample t-test

```r
#Two-sample t-test

#load libraries
library(tidyverse)
library(car)
library(moments)
library(ggplot2)

#test barnacle recruitment in open ocean vs. bays

#load data
barnacledata <- read.csv("Data/Barnacle data.csv")
barnacledata

#look at variance
#openocean
openvar <- barnacledata %>%
  filter(Site == "openocean") %>%
  summarize(var=var(recruitment, na.rm=TRUE))
openvar
#bays
bayvar <- barnacledata %>%
  filter(Site == "bay") %>%
  summarize(var=var(recruitment, na.rm=TRUE))
bayvar
#they're not the same, but within the same order of magnitude

#try a Bartlett test
bartlett.test(recruitment~Site, data=barnacledata) #P > 0.05

#check normality
#make a histogram
hist(barnacledata$recruitment, main= '') #very right skewed
shapiro.test(barnacledata$recruitment) #P << 0.05
skewness(barnacledata$recruitment) #1.01
kurtosis(barnacledata$recruitment) #4.10
qqPlot(barnacledata$recruitment)
#data is not normally distributed

#try a sqrt transformation (the data has zeros)
barnacledata$sqrtdata <- sqrt(barnacledata$recruitment)
barnacledata

#check normality
#make histogram
hist(barnacledata$sqrtdata, main = '') #very right skewed
shapiro.test(barnacledata$sqrtdata) #P = 0.003, this is better but still < 0.05
skewness(barnacledata$sqrtdata) #good
kurtosis(barnacledata$sqrtdata) #good
qqPlot(barnacledata$sqrtdata)
#data is still not normally distributed but the transformation improved it

#try a log transformation
#all the data get a +1 since there are zeros
barnacledata$plusone <- (barnacledata$recruitment)+1
barnacledata
barnacledata$logdata <- log(barnacledata$plusone)
barnacledata

#check normality
#make histogram
hist(barnacledata$logdata) #very weird
shapiro.test(barnacledata$logdata) #P << 0.05, this is worse
#not even gonna bother looking at skewness & kurtosis
qqPlot(barnacledata$logdata)
#log transformation doesn't help at all

#run two sample t-tests
#for a two sample test, the groups have to be specified as independetvar~categoricalvar
#with sqrt transformed data
t.test(sqrtdata~Site, var.equal=TRUE, data=barnacledata)
#with regular data
t.test(recruitment~Site, var.equal=TRUE, data=barnacledata)
#notice that these tests assume equal variance "var.equal = TRUE"
#if they weren't, instead run Welch's test by setting "var.equal = FALSE"
t.test(sqrtdata~Site, var.equal=FALSE, data=barnacledata)
#notice that the df and p-value are not the same

#try a Mann-Whitney U-test with raw data
wilcox.test(recruitment~Site, data=barnacledata)

#make a barplot with means +/- standard error bars
#first put together summary stats
graphdata <- barnacledata %>%
  group_by(Site) %>%
  summarize(meanrecruitment=mean(recruitment, na.rm=TRUE), se=sd(recruitment, na.rm=TRUE)/sqrt(30))
graphdata

#ggplot is best for making a plot with error bars

#first look at how to do this step-by-step
#notice that the object p (the plot) can be built upon with p+'s

#define what x and y variables to use to make the graph
p<-ggplot(graphdata, aes(x=Site, y=meanrecruitment)) 
#gts rid of the gray background and the useless gridlines
p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))
#specify what colors should the bars be. Why not think Dodger blue?
p<-p + geom_bar(stat="identity", color="dodgerblue2", fill="dodgerblue2", position="dodge", size=0.6) 
#label the axes
p<-p + labs(x="Site", y="Recruitment", size = 1) 
#add error bars
p<-p + geom_errorbar(aes(ymax=meanrecruitment+se, ymin=meanrecruitment-se), position=position_dodge(0.9), width=0.1) 
#print the final product
print(p)

#it's possible to do all of the above by instead adding the parts together step-wise
ggplot(graphdata, aes(x=Site, y=meanrecruitment)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_bar(stat="identity", color="dodgerblue2", fill="dodgerblue2", position="dodge", size=0.6) +
  labs(x="Site", y="Recruitment", size = 1) +
  geom_errorbar(aes(ymax=meanrecruitment+se, ymin=meanrecruitment-se), position=position_dodge(0.9), width=0.1) 

```

---

