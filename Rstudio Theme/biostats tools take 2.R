# Auto-sync background process launcher
.First <- function() {
  if (file.exists("auto_push.sh") && Sys.info()["sysname"] == "Windows") {
    # Full absolute path to PowerShell prevents the 'not found' error
    ps_path <- "C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe"
    
    if (file.exists(ps_path)) {
      cmd <- sprintf('"%s" -WindowStyle Hidden -Command "Start-Process bash -ArgumentList \\"./auto_push.sh\\" -WindowStyle Hidden"', ps_path)
      system(cmd, wait = FALSE, show.output.on.console = FALSE)
      message("Auto-sync watcher loaded in background.")
    } else {
      # Fallback directly through cmd.exe if PowerShell binary is elsewhere
      system('cmd.exe /c start /b bash ./auto_push.sh', wait = FALSE, show.output.on.console = FALSE)
      message("Auto-sync watcher loaded via cmd.")
    }
  }
}

# Manual sync helper
sync <- function(message = paste("Update:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))) {
  system2("git", c("add", "-A"))
  system2("git", c("commit", "-m", shQuote(message)))
  push_status <- system2("git", c("push", "origin", "HEAD"))
  if (push_status == 0) {
    message("Git push complete!")
  } else {
    warning("Git push encountered an issue.")
  }
}
--------------------------------------------------------------------------------------
  
  
  
  
  # ================================================================
# Joe's BIOL 330 BIOSTATISTICS TOOLKIT
# ================================================================
# This is the single master file to source from .Rprofile.
# Recommended startup command:
# source("C:/Rstudio/Biostatistics/biostats_toolkit.R")
# ================================================================


# ---------------------------------------------------------------
# 1. PACKAGE SETUP
# ---------------------------------------------------------------

load_biostats_packages <- function(install_missing = FALSE) {
  required_packages <- c("car")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  
  if (length(missing_packages) > 0 && install_missing) {
    install.packages(missing_packages)
    missing_packages <- missing_packages[
      !vapply(missing_packages, requireNamespace, logical(1), quietly = TRUE)
    ]
  }
  
  if (length(missing_packages) == 0) {
    message("Biostatistics toolkit loaded. Optional package car is available.")
  } else {
    message(
      "Biostatistics toolkit loaded. Optional package(s) unavailable: ",
      paste(missing_packages, collapse = ", "),
      ". Levene's tests will be skipped."
    )
  }
  
  invisible(missing_packages)
}


# ---------------------------------------------------------------
# 2. DATASET INSPECTION
# ---------------------------------------------------------------

inspect_df <- function(data,
                       digits = 3,
                       print_report = TRUE,
                       variable_types = NULL) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  
  if (!is.null(variable_types)) {
    if (is.null(names(variable_types))) {
      stop("variable_types must be a named character vector.")
    }
    unknown_types <- setdiff(variable_types, c("nominal", "ordinal"))
    if (length(unknown_types) > 0) {
      stop("variable_types may only contain 'nominal' or 'ordinal'.")
    }
  }
  
  classify_variable <- function(variable, x) {
    if (!is.null(variable_types) && variable %in% names(variable_types)) {
      return(c("categorical", variable_types[[variable]]))
    }
    if (is.ordered(x)) return(c("categorical", "ordinal"))
    if (is.factor(x) || is.character(x) || is.logical(x)) {
      return(c("categorical", "nominal"))
    }
    if (is.numeric(x)) return(c("quantitative", "numeric"))
    c("other", "other")
  }
  
  summary_table <- do.call(rbind, lapply(names(data), function(variable) {
    x <- data[[variable]]
    classification <- classify_variable(variable, x)
    
    if (is.numeric(x)) {
      finite_x <- x[is.finite(x)]
      n <- length(finite_x)
      mean_value <- if (n > 0) mean(finite_x) else NA_real_
      sd_value <- if (n > 1) sd(finite_x) else NA_real_
      standard_error_value <- if (n > 1) sd_value / sqrt(n) else NA_real_
      
      data.frame(
        variable = variable,
        category = classification[1],
        measurement_scale = classification[2],
        class = paste(class(x), collapse = "/"),
        n = n,
        missing = sum(is.na(x)),
        unique_values = length(unique(x[!is.na(x)])),
        minimum = if (n > 0) min(finite_x) else NA_real_,
        maximum = if (n > 0) max(finite_x) else NA_real_,
        mean = mean_value,
        median = if (n > 0) median(finite_x) else NA_real_,
        sd = sd_value,
        variance = if (n > 1) var(finite_x) else NA_real_,
        standard_error = standard_error_value,
        SEM = standard_error_value,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        variable = variable,
        category = classification[1],
        measurement_scale = classification[2],
        class = paste(class(x), collapse = "/"),
        n = sum(!is.na(x)),
        missing = sum(is.na(x)),
        unique_values = length(unique(x[!is.na(x)])),
        minimum = NA_real_, maximum = NA_real_, mean = NA_real_,
        median = NA_real_, sd = NA_real_, variance = NA_real_,
        standard_error = NA_real_, SEM = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }))
  
  numeric_names <- c(
    "minimum", "maximum", "mean", "median", "sd", "variance",
    "standard_error", "SEM"
  )
  summary_table[numeric_names] <- lapply(
    summary_table[numeric_names],
    function(x) round(x, digits)
  )
  
  if (print_report) {
    cat("\n====================================\n")
    cat("DATASET INSPECTION\n")
    cat("====================================\n")
    cat("Rows:", nrow(data), " | Columns:", ncol(data), "\n\n")
    
    cat("Categorical / nominal variables:\n")
    print(summary_table[summary_table$measurement_scale == "nominal",
                        c("variable", "class", "n", "missing", "unique_values")],
          row.names = FALSE)
    
    cat("\nCategorical / ordinal variables:\n")
    print(summary_table[summary_table$measurement_scale == "ordinal",
                        c("variable", "class", "n", "missing", "unique_values")],
          row.names = FALSE)
    
    cat("\nQuantitative / numeric variables:\n")
    print(summary_table[summary_table$category == "quantitative", ],
          row.names = FALSE)
    
    cat("\nNote: standard_error and SEM are the same here. SEM means standard error of the mean.\n")
  }
  
  invisible(summary_table)
}


# ---------------------------------------------------------------
# 3. TRANSFORMATIONS
# ---------------------------------------------------------------

add_transformations <- function(data, variable, save_to = NULL) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!variable %in% names(data)) stop("The variable was not found.")
  if (!is.numeric(data[[variable]])) stop("The variable must be numeric.")
  
  x <- data[[variable]]
  
  sqrt_values <- rep(NA_real_, length(x))
  nonnegative <- !is.na(x) & x >= 0
  sqrt_values[nonnegative] <- sqrt(x[nonnegative])
  
  log_values <- rep(NA_real_, length(x))
  positive <- !is.na(x) & x > 0
  log_values[positive] <- log(x[positive])
  
  log10_values <- rep(NA_real_, length(x))
  log10_values[positive] <- log10(x[positive])
  
  data[[paste0(variable, "_sqrt")]] <- sqrt_values
  data[[paste0(variable, "_square")]] <- x^2
  data[[paste0(variable, "_log")]] <- log_values
  data[[paste0(variable, "_log10")]] <- log10_values
  
  if (!is.null(save_to)) {
    write.csv(data, save_to, row.names = FALSE)
    cat("Updated data saved to:", save_to, "\n")
  }
  
  data
}


# ---------------------------------------------------------------
# 4. Q-Q PLOTS WITH CONFIDENCE ENVELOPES
# ---------------------------------------------------------------

qq_plot_with_ci <- function(x,
                            conf = 0.95,
                            simulations = 1000,
                            main = "Q-Q Plot with Confidence Envelope") {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 3) stop("At least 3 finite observations are needed.")
  if (sd(x) == 0) stop("Q-Q plot cannot be created when all values are identical.")
  
  theoretical_values <- qnorm(ppoints(n))
  ordered_values <- sort(x)
  
  simulated_values <- matrix(
    rnorm(n * simulations, mean = mean(x), sd = sd(x)),
    nrow = simulations,
    ncol = n
  )
  simulated_values <- t(apply(simulated_values, 1, sort))
  
  lower_probability <- (1 - conf) / 2
  upper_probability <- 1 - lower_probability
  lower_limit <- apply(simulated_values, 2, quantile,
                       probs = lower_probability)
  upper_limit <- apply(simulated_values, 2, quantile,
                       probs = upper_probability)
  
  plot(theoretical_values, ordered_values,
       main = main,
       xlab = "Theoretical normal quantiles",
       ylab = "Observed quantiles",
       pch = 19)
  lines(theoretical_values, lower_limit,
        col = "blue", lty = 2, lwd = 2)
  lines(theoretical_values, upper_limit,
        col = "blue", lty = 2, lwd = 2)
  qqline(ordered_values, col = "red", lwd = 2)
  
  legend("topleft",
         legend = c("Reference line", paste0(conf * 100, "% envelope")),
         col = c("red", "blue"),
         lty = c(1, 2),
         lwd = 2,
         bty = "n")
  
  invisible(list(
    theoretical = theoretical_values,
    observed = ordered_values,
    lower = lower_limit,
    upper = upper_limit
  ))
}


# ---------------------------------------------------------------
# 5. GENERAL NORMALITY CHECK
# ---------------------------------------------------------------

assumption_report <- function(data,
                              response = NULL,
                              group = NULL,
                              alpha = 0.05,
                              plots = TRUE) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  
  if (is.null(response)) {
    response <- names(data)[vapply(data, is.numeric, logical(1))]
  }
  if (length(response) == 0) stop("No numeric response variables were found.")
  
  if (!is.null(group) && !group %in% names(data)) {
    stop("The grouping variable was not found.")
  }
  
  results <- list()
  
  for (variable in response) {
    if (!variable %in% names(data)) {
      stop("Response variable was not found: ", variable)
    }
    
    selected_columns <- c(variable, group)
    working_data <- data[, selected_columns, drop = FALSE]
    working_data <- working_data[complete.cases(working_data), , drop = FALSE]
    values <- working_data[[variable]]
    
    cat("\n====================================\n")
    cat("Response variable:", variable, "\n")
    cat("Complete observations:", length(values), "\n")
    cat("====================================\n")
    
    overall_shapiro <- NULL
    if (length(values) >= 3 && length(values) <= 5000) {
      overall_shapiro <- shapiro.test(values)
      cat("\nOverall Shapiro-Wilk test:\n")
      print(overall_shapiro)
    }
    
    group_shapiro <- NULL
    levene_result <- NULL
    
    if (!is.null(group)) {
      group_shapiro <- lapply(split(values, working_data[[group]]), function(x) {
        if (length(x) >= 3 && length(x) <= 5000) shapiro.test(x) else NULL
      })
      
      cat("\nShapiro-Wilk tests by group:\n")
      for (group_name in names(group_shapiro)) {
        cat("\nGroup:", group_name, "\n")
        if (is.null(group_shapiro[[group_name]])) {
          cat("Unavailable: group needs 3 to 5000 observations.\n")
        } else {
          print(group_shapiro[[group_name]])
        }
      }
      
      if (requireNamespace("car", quietly = TRUE)) {
        levene_data <- data.frame(
          response = values,
          group = as.factor(working_data[[group]])
        )
        levene_result <- car::leveneTest(response ~ group, data = levene_data)
        cat("\nLevene's test for equal variances:\n")
        print(levene_result)
      }
    }
    
    if (plots && length(values) >= 3) {
      qq_plot_with_ci(values, main = paste("Q-Q Plot:", variable))
    }
    
    results[[variable]] <- list(
      n = length(values),
      overall_shapiro = overall_shapiro,
      group_shapiro = group_shapiro,
      levene = levene_result
    )
  }
  
  invisible(results)
}


# ---------------------------------------------------------------
# 6. TRANSFORMATION WORKFLOW
# ---------------------------------------------------------------

transformation_workflow <- function(data,
                                    variable,
                                    alpha = 0.05,
                                    plots = TRUE,
                                    save_to = NULL) {
  transformed_data <- add_transformations(data, variable, save_to = save_to)
  
  variables_to_check <- c(
    variable,
    paste0(variable, "_sqrt"),
    paste0(variable, "_square"),
    paste0(variable, "_log"),
    paste0(variable, "_log10")
  )
  
  test_results <- lapply(variables_to_check, function(current_variable) {
    values <- transformed_data[[current_variable]]
    values <- values[is.finite(values)]
    
    if (length(values) < 3 || length(values) > 5000) {
      return(data.frame(
        variable = current_variable,
        n = length(values),
        shapiro_p_value = NA_real_,
        interpretation = "Shapiro-Wilk unavailable",
        stringsAsFactors = FALSE
      ))
    }
    
    shapiro_result <- shapiro.test(values)
    interpretation <- if (shapiro_result$p.value < alpha) {
      "Evidence against normality"
    } else {
      "No significant evidence against normality"
    }
    
    if (plots) {
      qq_plot_with_ci(values, main = paste("Q-Q Plot:", current_variable))
    }
    
    data.frame(
      variable = current_variable,
      n = length(values),
      shapiro_p_value = shapiro_result$p.value,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  })
  
  test_results <- do.call(rbind, test_results)
  cat("\nTransformation normality comparison:\n")
  print(test_results, row.names = FALSE)
  
  invisible(list(
    data = transformed_data,
    normality_results = test_results
  ))
}

compare_transformations <- function(data,
                                    variable,
                                    alpha = 0.05,
                                    plots = TRUE,
                                    save_to = NULL) {
  transformation_workflow(
    data = data,
    variable = variable,
    alpha = alpha,
    plots = plots,
    save_to = save_to
  )
}


# ---------------------------------------------------------------
# 7. ANOVA WORKFLOW
# ---------------------------------------------------------------

anova_workflow <- function(data,
                           formula,
                           alpha = 0.05,
                           plots = TRUE,
                           run_tukey = TRUE) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!inherits(formula, "formula")) stop("formula must look like y ~ group.")
  
  model <- lm(formula, data = data, na.action = na.omit)
  model_frame <- model.frame(model)
  predictor_names <- all.vars(delete.response(terms(model)))
  response_name <- all.vars(formula)[1]
  residual_values <- residuals(model)
  fitted_values <- fitted(model)
  
  if (length(predictor_names) == 0) {
    stop("The formula must include at least one predictor.")
  }
  
  numeric_predictors <- predictor_names[
    vapply(model_frame[predictor_names], is.numeric, logical(1))
  ]
  if (length(numeric_predictors) > 0) {
    warning(
      "These predictors are numeric, so this may be regression rather than ANOVA: ",
      paste(numeric_predictors, collapse = ", ")
    )
  }
  
  cat("\n====================================\n")
  cat("ANOVA WORKFLOW\n")
  cat("====================================\n")
  cat("Model: ", deparse(formula), "\n", sep = "")
  cat("Complete observations used:", nrow(model_frame), "\n")
  
  anova_table <- anova(model)
  cat("\nANOVA table:\n")
  print(anova_table)
  
  cat("\nANOVA interpretation:\n")
  anova_p_values <- anova_table[["Pr(>F)"]]
  for (term_number in seq_along(anova_p_values)) {
    term_name <- rownames(anova_table)[term_number]
    p_value <- anova_p_values[term_number]
    
    if (is.na(p_value)) {
      next
    }
    
    if (p_value < alpha) {
      cat(
        "-", term_name, ": p =", format.pval(p_value),
        "-> evidence that means differ for this term.\n"
      )
    } else {
      cat(
        "-", term_name, ": p =", format.pval(p_value),
        "-> no significant evidence that means differ for this term.\n"
      )
    }
  }
  
  residual_diagnostics <- data.frame(
    observation = seq_along(residual_values),
    fitted_value = fitted_values,
    residual = residual_values,
    standardized_residual = rstandard(model)
  )
  
  cat("\nResidual diagnostics:\n")
  cat("Residuals are observed values minus model-predicted values.\n")
  cat("A residual near 0 means the model predicted that observation well.\n")
  cat("Large positive or negative residuals are observations the model predicted poorly.\n")
  cat("\nFirst residuals:\n")
  print(head(residual_diagnostics, 10), row.names = FALSE)
  cat("\nResidual summary:\n")
  print(summary(residual_values))
  
  shapiro_result <- NULL
  if (length(residual_values) >= 3 && length(residual_values) <= 5000) {
    shapiro_result <- shapiro.test(residual_values)
    cat("\nShapiro-Wilk test on model residuals:\n")
    print(shapiro_result)
    
    if (shapiro_result$p.value < alpha) {
      cat("Interpretation: evidence against normally distributed residuals.\n")
      cat("This suggests inspecting the Q-Q plot and considering transformation or a different model.\n")
    } else {
      cat("Interpretation: no significant evidence against normal residuals.\n")
      cat("This assumption appears reasonable based on this test.\n")
    }
  }
  
  grouping_data <- model_frame[predictor_names]
  group_label <- do.call(
    interaction,
    c(grouping_data, list(drop = TRUE, sep = " × "))
  )
  
  response_values <- model_frame[[response_name]]
  levene_result <- NULL
  if (requireNamespace("car", quietly = TRUE)) {
    levene_data <- data.frame(response = response_values, group = group_label)
    levene_result <- car::leveneTest(response ~ group, data = levene_data)
    cat("\nLevene's test for equal variances:\n")
    print(levene_result)
    levene_p_value <- levene_result[["Pr(>F)"]][1]
    if (levene_p_value < alpha) {
      cat("Interpretation: evidence that group variances are unequal.\n")
    } else {
      cat("Interpretation: no significant evidence that group variances differ.\n")
    }
  } else {
    cat("\nLevene's test skipped. To enable it, run:\n")
    cat('install.packages("car")\n')
  }
  
  split_values <- split(response_values, group_label)
  group_summary <- do.call(rbind, lapply(names(split_values), function(group_name) {
    x <- split_values[[group_name]]
    n <- length(x)
    mean_value <- mean(x)
    sd_value <- if (n > 1) sd(x) else NA_real_
    se_value <- if (n > 1) sd_value / sqrt(n) else NA_real_
    margin <- if (n > 1) qt(.975, df = n - 1) * se_value else NA_real_
    
    data.frame(
      group = group_name,
      n = n,
      mean = mean_value,
      sd = sd_value,
      se = se_value,
      lower_95_CI = mean_value - margin,
      upper_95_CI = mean_value + margin
    )
  }))
  rownames(group_summary) <- NULL
  
  cat("\nGroup summary:\n")
  print(group_summary, row.names = FALSE)
  
  anova_ss <- anova_table[["Sum Sq"]]
  eta_squared <- sum(anova_ss[-length(anova_ss)]) / sum(anova_ss)
  cat("\nOverall eta-squared:", round(eta_squared, 4), "\n")
  eta_interpretation <- if (eta_squared < 0.01) {
    "very small effect"
  } else if (eta_squared < 0.06) {
    "small effect"
  } else if (eta_squared < 0.14) {
    "medium effect"
  } else {
    "large effect"
  }
  cat("Approximate effect-size interpretation:", eta_interpretation, "\n")
  
  tukey_result <- NULL
  if (run_tukey && length(predictor_names) == 1 &&
      is.factor(model_frame[[predictor_names]])) {
    tukey_result <- TukeyHSD(aov(model))
    cat("\nTukey HSD post-hoc comparisons:\n")
    print(tukey_result)
    cat("Interpretation: Tukey compares specific pairs of groups after the overall ANOVA.\n")
    cat("Pairs with adjusted p-values below", alpha, "show evidence of differing means.\n")
  } else if (run_tukey) {
    cat("\nTukey HSD was not automatically run because this is not a simple one-factor ANOVA.\n")
  }
  
  if (plots) {
    old_par <- par(mfrow = c(1, 2))
    on.exit(par(old_par), add = TRUE)
    
    qq_plot_with_ci(
      residual_values,
      main = "Q-Q Plot of ANOVA Residuals"
    )
    plot(fitted_values, residual_values,
         main = "Residuals vs. Fitted Values",
         xlab = "Fitted values",
         ylab = "Residuals")
    abline(h = 0, col = "red", lwd = 2)
    
    cat("\nPlot interpretation:\n")
    cat("- Q-Q plot: residuals should generally follow the red line and remain near the blue envelope.\n")
    cat("- Residuals-versus-fitted plot: residuals should be randomly scattered around zero without a funnel or curve.\n")
  }
  
  invisible(list(
    model = model,
    anova = anova_table,
    residual_diagnostics = residual_diagnostics,
    shapiro = shapiro_result,
    levene = levene_result,
    group_summary = group_summary,
    eta_squared = eta_squared,
    tukey = tukey_result
  ))
}


anova_assumption_report <- function(data,
                                    formula,
                                    alpha = 0.05,
                                    plots = TRUE) {
  anova_workflow(
    data = data,
    formula = formula,
    alpha = alpha,
    plots = plots,
    run_tukey = FALSE
  )
}


# Load optional packages when this master file is sourced.
load_biostats_packages()

