# Personal biostatistics toolkit

add_transformations <- function(data, variable) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!variable %in% names(data)) stop("The variable was not found.")
  if (!is.numeric(data[[variable]])) stop("The variable must be numeric.")

  x <- data[[variable]]
  data[[paste0(variable, "_sqrt")]] <- ifelse(x >= 0, sqrt(x), NA_real_)
  data[[paste0(variable, "_square")]] <- x^2
  data[[paste0(variable, "_log")]] <- ifelse(x > 0, log(x), NA_real_)
  data[[paste0(variable, "_log10")]] <- ifelse(x > 0, log10(x), NA_real_)

  data
}

anova_workflow <- function(data,
                           formula,
                           alpha = 0.05,
                           plots = TRUE,
                           run_tukey = TRUE) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!inherits(formula, "formula")) stop("formula must look like y ~ group.")

  model <- lm(formula, data = data, na.action = na.omit)
  model_frame <- model.frame(model)
  model_terms <- attr(terms(model), "term.labels")
  predictor_names <- all.vars(delete.response(terms(model)))
  response_name <- all.vars(formula)[1]
  residual_values <- residuals(model)
  fitted_values <- fitted(model)

  if (length(predictor_names) == 0) {
    stop("The formula must include at least one predictor.")
  }

  cat("\n====================================\n")
  cat("ANOVA WORKFLOW\n")
  cat("====================================\n")
  cat("Model: ", deparse(formula), "\n", sep = "")
  cat("Response variable: ", response_name, "\n", sep = "")
  cat("Complete observations used: ", nrow(model_frame), "\n", sep = "")

  # ANOVA table
  anova_table <- anova(model)
  cat("\nANOVA table:\n")
  print(anova_table)

  # Residual normality
  shapiro_result <- NULL
  if (length(residual_values) >= 3 && length(residual_values) <= 5000) {
    shapiro_result <- shapiro.test(residual_values)
    cat("\nShapiro-Wilk test on model residuals:\n")
    print(shapiro_result)
    if (shapiro_result$p.value < alpha) {
      cat("Interpretation: evidence against normal residuals.\n")
    } else {
      cat("Interpretation: no significant evidence against normal residuals.\n")
    }
  } else {
    cat("\nShapiro-Wilk test skipped: sample size must be between 3 and 5000.\n")
  }

  # Equal-variance test
  levene_result <- NULL
  if (requireNamespace("car", quietly = TRUE)) {
    levene_result <- car::leveneTest(formula, data = model_frame)
    cat("\nLevene's test for equal variances:\n")
    print(levene_result)
  } else {
    cat("\nLevene's test skipped. Install the car package with:\n")
    cat('install.packages("car")\n')
  }

  # Group summary using all predictor combinations
  grouping_data <- model_frame[predictor_names]
  group_label <- do.call(
    interaction,
    c(grouping_data, list(drop = TRUE, sep = " × "))
  )
  summary_data <- data.frame(
    group = group_label,
    response = model_frame[[response_name]]
  )

  group_summary <- do.call(rbind, lapply(split(summary_data$response,
                                                summary_data$group), function(x) {
    n <- length(x)
    mean_value <- mean(x)
    sd_value <- sd(x)
    data.frame(
      n = n,
      mean = mean_value,
      sd = sd_value,
      se = sd_value / sqrt(n),
      lower_95_CI = mean_value - qt(.975, df = n - 1) * sd_value / sqrt(n),
      upper_95_CI = mean_value + qt(.975, df = n - 1) * sd_value / sqrt(n)
    )
  }))
  group_summary$group <- rownames(group_summary)
  rownames(group_summary) <- NULL
  group_summary <- group_summary[, c("group", "n", "mean", "sd", "se",
                                     "lower_95_CI", "upper_95_CI")]

  cat("\nGroup summary:\n")
  print(group_summary, row.names = FALSE)

  # Effect size: eta squared for the overall model
  anova_ss <- anova_table[["Sum Sq"]]
  eta_squared <- sum(anova_ss[-length(anova_ss)]) / sum(anova_ss)
  cat("\nOverall eta-squared:", round(eta_squared, 4), "\n")

  # Tukey post-hoc test for a simple one-factor ANOVA
  tukey_result <- NULL
  if (run_tukey && length(predictor_names) == 1 &&
      is.factor(model_frame[[predictor_names]])) {
    tukey_result <- TukeyHSD(aov(model))
    cat("\nTukey HSD post-hoc comparisons:\n")
    print(tukey_result)
  } else if (run_tukey) {
    cat("\nTukey HSD was not automatically run because this is not a simple one-factor ANOVA.\n")
  }

  # Diagnostic plots
  if (plots) {
    old_par <- par(mfrow = c(1, 2))
    on.exit(par(old_par), add = TRUE)

    qqnorm(residual_values, main = "Q-Q Plot of ANOVA Residuals")
    qqline(residual_values, col = "red", lwd = 2)

    plot(fitted_values, residual_values,
         main = "Residuals vs. Fitted Values",
         xlab = "Fitted values", ylab = "Residuals")
    abline(h = 0, col = "red", lwd = 2)
  }

  invisible(list(
    model = model,
    anova = anova_table,
    shapiro = shapiro_result,
    levene = levene_result,
    group_summary = group_summary,
    eta_squared = eta_squared,
    tukey = tukey_result
  ))
}
