git_bin <- "C:/Program Files/Git/bin/git.exe"
if (!file.exists(git_bin)) {
  git_bin <- Sys.which("git")
}

# Set author environment variables directly for this process
Sys.setenv(GIT_AUTHOR_NAME = "Joseph Goldenberg")
Sys.setenv(GIT_AUTHOR_EMAIL = "ijosephgold@gmail.com")
Sys.setenv(GIT_COMMITTER_NAME = "Joseph Goldenberg")
Sys.setenv(GIT_COMMITTER_EMAIL = "ijosephgold@gmail.com")

repeat {
  status <- system2(git_bin, c("status", "--porcelain"), stdout = TRUE)
  if (length(status) > 0) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    system2(git_bin, c("add", "-A"))
    system2(git_bin, c("commit", "-m", shQuote(paste("Auto-sync:", timestamp))))
    system2(git_bin, c("push", "origin", "master"))
  }
  Sys.sleep(5)
}
