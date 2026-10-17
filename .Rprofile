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

# Dataset inspection helper
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
