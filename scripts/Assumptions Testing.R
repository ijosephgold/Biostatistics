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
    "Could not test: group size must be 3–5000",
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