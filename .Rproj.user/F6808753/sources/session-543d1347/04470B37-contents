choose_column <- function(data, candidates, question) {
  if (length(candidates) == 0) {
    stop("No suitable columns were found for this question.")
  }
  
  choices <- vapply(candidates, function(name) {
    column_number <- match(name, names(data))
    values <- unique(na.omit(data[[name]]))
    examples <- paste(head(format(values, trim = TRUE), 3), collapse = ", ")
    
    paste0(
      "Data column ", column_number, " — ", name,
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
    "Could not test: sample size must be 3–5000",
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