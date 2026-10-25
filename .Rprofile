# Biostatistics project startup
# Paths are relative to this repository so the project works on both laptops.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

toolkit_file <- file.path(project_root, "biostats_toolkit.R")
if (file.exists(toolkit_file)) {
  source(toolkit_file, local = .GlobalEnv)
} else {
  message("Biostatistics toolkit file was not found: ", toolkit_file)
}

sync_file <- file.path(project_root, "run_sync.R")
if (file.exists(sync_file)) {
  source(sync_file, local = .GlobalEnv)
} else {
  message("Automatic sync script was not found: ", sync_file)
}
