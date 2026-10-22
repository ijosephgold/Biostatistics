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