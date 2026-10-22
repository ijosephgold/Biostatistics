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

