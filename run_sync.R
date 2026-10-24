# ============================================================
# BIOL 330 AUTOMATIC GITHUB SYNC
# ============================================================

project_dir <- "C:/Rstudio/Biostatistics"

git_bin <- "C:/Program Files/Git/bin/git.exe"

if (!file.exists(git_bin)) {
  git_bin <- Sys.which("git")
}

if (git_bin == "" || !file.exists(git_bin)) {
  stop("Git could not be found.")
}

if (!dir.exists(project_dir)) {
  stop("The Biostatistics project folder could not be found.")
}

# Confirm that this folder is actually a Git repository
repo_check <- system2(
  git_bin,
  c("-C", project_dir, "rev-parse", "--is-inside-work-tree"),
  stdout = TRUE,
  stderr = FALSE
)

if (!any(trimws(repo_check) == "true")) {
  stop(
    "This folder is not connected to Git. ",
    "A .git folder or GitHub connection is missing."
  )
}

# Git author information
Sys.setenv(
  GIT_AUTHOR_NAME = "Joseph Goldenberg",
  GIT_AUTHOR_EMAIL = "ijosephgold@gmail.com",
  GIT_COMMITTER_NAME = "Joseph Goldenberg",
  GIT_COMMITTER_EMAIL = "ijosephgold@gmail.com"
)

cat("Automatic GitHub sync started.\n")
cat("Watching:", project_dir, "\n")

repeat {
  
  status <- system2(
    git_bin,
    c("-C", project_dir, "status", "--porcelain"),
    stdout = TRUE,
    stderr = FALSE
  )
  
  # Only commit and push when actual changes exist
  if (length(status) > 0 && any(nzchar(status))) {
    
    timestamp <- format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
    
    system2(
      git_bin,
      c("-C", project_dir, "add", "-A")
    )
    
    system2(
      git_bin,
      c(
        "-C", project_dir,
        "commit",
        "-m",
        shQuote(paste("Automatic sync:", timestamp))
      )
    )
    
    push_result <- system2(
      git_bin,
      c("-C", project_dir, "push", "origin", "master"),
      stdout = TRUE,
      stderr = TRUE
    )
    
    cat("Synced:", timestamp, "\n")
    print(push_result)
  }
  
  # Check once per minute
  Sys.sleep(60)
}