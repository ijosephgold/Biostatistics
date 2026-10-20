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
        "Data column ", column_number, " — ", name,
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
