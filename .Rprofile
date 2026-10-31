# ============================================================
# AUTOMATIC BIOSTATISTICS TOOLKIT LOADER
# ============================================================

# Load the toolkit
toolkit_file <- "C:/Rstudio/Biostatistics/biostats_toolkit.R"

if (file.exists(toolkit_file)) {
  source(toolkit_file, local = .GlobalEnv)
}

# Load the biostatistics packages
if (exists("load_biostats_packages", mode = "function")) {
  load_biostats_packages(install_missing = FALSE)
}
