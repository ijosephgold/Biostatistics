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
  cat("========================================================================================\n")
  cat("                                    DATASET OVERVIEW                                    \n")
  cat("========================================================================================\n")
  cat(sprintf("Dimensions  : %d rows x %d columns\n", nrow(df), ncol(df)))
  total_cells <- prod(dim(df))
  total_na <- sum(is.na(df))
  cat(sprintf("Total NA    : %d (%.1f%% of all cells)\n", total_na, (total_na / total_cells) * 100))
  cat(sprintf("Completeness: %d complete rows (%.1f%%)\n\n", 
              sum(complete.cases(df)), 
              (sum(complete.cases(df)) / nrow(df)) * 100))
  
  calc_mode <- function(x) {
    clean_x <- x[!is.na(x)]
    if (length(clean_x) == 0) return(NA)
    ux <- unique(clean_x)
    tab <- tabulate(match(clean_x, ux))
    ux[tab == max(tab)][1]
  }
  
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) > 0) {
    cat("--- Numeric Variables ---\n")
    num_list <- lapply(num_cols, function(col) {
      vals <- df[[col]]
      clean <- vals[!is.na(vals)]
      n_clean <- length(clean)
      
      data.frame(
        Variable = col,
        N = n_clean,
        Unique = length(unique(clean)),
        NonZero_N = sum(clean != 0),
        Mean = round(mean(clean), digits),
        SD = round(sd(clean), digits),
        SEM = round(sd(clean) / sqrt(n_clean), digits),
        Var = round(var(clean), digits),
        Median = round(median(clean), digits),
        IQR = round(IQR(clean), digits),
        Min = round(min(clean), digits),
        Max = round(max(clean), digits),
        Mode = round(calc_mode(clean), digits),
        stringsAsFactors = FALSE
      )
    })
    print(do.call(rbind, num_list), row.names = FALSE)
    cat("\n")
  }
  
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
  cat("========================================================================================\n")
}
