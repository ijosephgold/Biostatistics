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
inspect_df <- function(df) {
  cat("--- Dataset Overview ---\n")
  cat("Dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
  cat("Total NA values:", sum(is.na(df)), "\n\n")
  cat("--- Column Breakdown ---\n")
  summary_info <- data.frame(
    Column = names(df),
    Type = sapply(df, function(x) class(x)[1]),
    Missing = sapply(df, function(x) sum(is.na(x))),
    Missing_Pct = paste0(round(sapply(df, function(x) mean(is.na(x))) * 100, 1), "%"),
    Unique = sapply(df, function(x) length(unique(x))),
    row.names = NULL
  )
  print(summary_info)
}
