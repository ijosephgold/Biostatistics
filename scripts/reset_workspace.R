reset_workspace <- function() {
  rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  gc(verbose = FALSE)
  cat("\014") # Clears the console screen
  message("Global environment cleared and memory freed.")
}