# ============================================================
# AUTOMATIC BIOSTATISTICS TOOLKIT LOADER
# ============================================================

# Load the toolkit
toolkit_file <- "C:/Rstudio/Biostatistics/biostats_toolkit.R"

if (file.exists(toolkit_file)) {
  source(toolkit_file, local = .GlobalEnv)
}

packages_to_attach <- c(
  "car",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "rstatix",
  "moments",
  "effectsize",
  "emmeans",
  "performance",
  "moments",
  "lme4",
  "tidyverse",
  "agricolae",
  "vegan",
  "rstatix",
  "afex",
  "ggpubr",
  "emmeans"
)

invisible(lapply(packages_to_attach, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  } else {
    message("Not installed yet; skipping: ", pkg)
  }
}))

rm(packages_to_attach)
# ==============================================================================
# AUTO-LOAD Gbiostats_toolkit.R FROM C:/Rstudio/Biostatistics
# ==============================================================================
local({
  # 1. Check C:/Rstudio/Biostatistics, current project folder, and ~
  candidates <- c(
    "C:/Rstudio/Biostatistics/Gbiostats_toolkit.R",
    "C:/Rstudio/Biostatistics/Gbiostats_toolkit.r",
    "Gbiostats_toolkit.R",
    "Gbiostats_toolkit.r",
    path.expand("~/Gbiostats_toolkit.R"),
    path.expand("~/Gbiostats_toolkit.r")
  )
  existing <- candidates[file.exists(candidates)]
  
  if (length(existing) > 0) {
    toolkit_path <- normalizePath(existing[1], winslash = "/")
    
    # 2. Prevent duplicate attachment if .Rprofile is re-sourced
    env_name <- "tools:Gbiostats"
    if (env_name %in% search()) {
      detach(env_name, character.only = TRUE)
    }
    
    # 3. Source into a dedicated environment and attach to search path
    biostat_env <- new.env(parent = globalenv())
    sys.source(toolkit_path, envir = biostat_env)
    attach(biostat_env, name = env_name, warn.conflicts = FALSE)
  } else {
    warning("Could not find Gbiostats_toolkit.R in C:/Rstudio/Biostatistics/")
  }
})
