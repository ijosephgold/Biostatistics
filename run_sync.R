# Automatic Biostatistics repository synchronization
#
# Sourcing this file starts a detached watcher and returns immediately.
# The watcher commits local coursework, fetches, integrates safely, then pushes.

biostats_sync_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script <- if (length(file_arg)) {
    sub("^--file=", "", file_arg[[1]])
  } else {
    file.path(getwd(), "run_sync.R")
  }
  normalizePath(dirname(script), winslash = "/", mustWork = TRUE)
}

biostats_sync_log <- function(root) file.path(root, ".git", "automatic-sync.log")
biostats_sync_lock <- function(root) file.path(root, ".git", "automatic-sync.lock")
biostats_sync_stop <- function(root) file.path(root, ".git", "automatic-sync.stop")

biostats_sync_write <- function(root, ...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ",
                 paste(..., collapse = " "))
  cat(line, "\n", file = biostats_sync_log(root), append = TRUE)
}

biostats_sync_git <- function(root, args, timeout = 60) {
  git <- Sys.which("git")
  if (!nzchar(git)) {
    candidates <- c(
      file.path(Sys.getenv("ProgramFiles"), "Git", "bin", "git.exe"),
      "C:/Program Files/Git/bin/git.exe"
    )
    candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
    if (length(candidates)) git <- candidates[[1]]
  }
  if (!nzchar(git)) {
    return(list(status = 127L,
                output = "Git was not found on PATH or in the standard Windows installation folder."))
  }

  quoted <- vapply(args, shQuote, character(1),
                   type = if (.Platform$OS.type == "windows") "cmd" else "sh")
  output <- tryCatch(
    suppressWarnings(system2(
      git, c("-C", shQuote(root,
        type = if (.Platform$OS.type == "windows") "cmd" else "sh"), quoted),
      stdout = TRUE, stderr = TRUE, stdin = "", timeout = timeout
    )),
    error = function(e) structure(conditionMessage(e), status = 127L)
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), output = as.character(output))
}

biostats_sync_report_git <- function(root, label, result) {
  if (length(result$output) && any(nzchar(result$output))) {
    for (line in result$output[nzchar(result$output)]) {
      biostats_sync_write(root, label, "|", line)
    }
  }
  if (result$status != 0L) {
    biostats_sync_write(root, label, "FAILED with exit code", result$status)
  }
  invisible(result$status == 0L)
}

biostats_sync_worker_running <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (is.na(pid) || pid < 1L) return(FALSE)
  tryCatch(isTRUE(tools::pskill(pid, signal = 0)),
           error = function(e) FALSE)
}

biostats_sync_start <- function(root = biostats_sync_root()) {
  if (!dir.exists(file.path(root, ".git"))) {
    # Worktrees use a .git file, so allow either form.
    if (!file.exists(file.path(root, ".git"))) {
      message("Automatic sync not started: no .git directory or file at ", root)
      return(invisible(FALSE))
    }
  }

  if (file.exists(biostats_sync_stop(root))) {
    message("Automatic sync is paused. Remove .git/automatic-sync.stop and reopen the project to resume.")
    return(invisible(FALSE))
  }

  lock <- biostats_sync_lock(root)
  if (dir.exists(lock)) {
    pid_file <- file.path(lock, "pid")
    old_pid <- if (file.exists(pid_file)) readLines(pid_file, warn = FALSE)[1] else NA_character_
    if (length(old_pid) && biostats_sync_worker_running(old_pid)) {
      message("Automatic GitHub sync is already running in the background.")
      return(invisible(TRUE))
    }
    unlink(lock, recursive = TRUE, force = TRUE)
  }

  dir.create(dirname(biostats_sync_log(root)), recursive = TRUE, showWarnings = FALSE)
  rscript <- file.path(R.home("bin"),
                       if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (!file.exists(rscript)) rscript <- Sys.which("Rscript")
  if (!nzchar(rscript) || !file.exists(rscript)) {
    message("Automatic sync not started: Rscript could not be found.")
    return(invisible(FALSE))
  }

  script <- file.path(root, "run_sync.R")
  log_file <- biostats_sync_log(root)
  worker_args <- c("--vanilla", shQuote(script,
    type = if (.Platform$OS.type == "windows") "cmd" else "sh"), "--worker")
  launch <- tryCatch(
    suppressWarnings(system2(
      rscript, worker_args, wait = FALSE, stdout = log_file,
      stderr = log_file, invisible = TRUE
    )),
    error = function(e) e
  )
  if (inherits(launch, "error")) {
    biostats_sync_write(root, "ERROR: could not start background watcher:",
                        conditionMessage(launch))
    message("Automatic sync could not start. See ", log_file)
    return(invisible(FALSE))
  }

  message("Automatic GitHub sync is running in the background. Log: ", log_file)
  invisible(TRUE)
}

biostats_sync_once <- function(root) {
  git <- function(...) biostats_sync_git(root, c(...))
  branch <- git("symbolic-ref", "--short", "HEAD")
  if (!biostats_sync_report_git(root, "Check current branch", branch)) return(FALSE)
  if (length(branch$output) != 1L || branch$output[[1]] != "master") {
    biostats_sync_write(root, "PAUSED: current branch is not master.")
    return(FALSE)
  }

  remote <- git("remote", "get-url", "origin")
  if (!biostats_sync_report_git(root, "Check origin", remote)) return(FALSE)

  in_progress <- git("rev-parse", "--git-path", "rebase-merge")
  rebase_path <- if (length(in_progress$output)) in_progress$output[[1]] else ""
  if (nzchar(rebase_path) && dir.exists(file.path(root, rebase_path))) {
    biostats_sync_write(root, "PAUSED: a rebase is already in progress; resolve it manually.")
    return(FALSE)
  }

  status <- git("status", "--porcelain", "--untracked-files=normal")
  if (!biostats_sync_report_git(root, "Read working-tree status", status)) return(FALSE)

  if (length(status$output) && any(nzchar(status$output))) {
    # Remove already-staged session edits from the index without changing their
    # working copies. Keep staged deletions so the documented untracking command
    # can be committed, but never include changing session-state contents.
    session_staged <- git("diff", "--cached", "--name-status", "--",
                          ".Rproj.user", ".RData", ".Rhistory", ".Ruserdata")
    if (session_staged$status == 0L && length(session_staged$output)) {
      for (line in session_staged$output[nzchar(session_staged$output)]) {
        fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
        if (length(fields) >= 2L && fields[[1]] != "D") {
          unstage <- git("reset", "-q", "--", fields[-1])
          if (!biostats_sync_report_git(root, "Unstage RStudio session state", unstage)) {
            return(FALSE)
          }
        }
      }
    } else if (session_staged$status != 0L) {
      biostats_sync_report_git(root, "Check staged RStudio session state", session_staged)
      return(FALSE)
    }

    # RStudio session/workspace state is excluded even if it was tracked in the past.
    add <- git("add", "-A", "--", ".",
               ":(exclude).Rproj.user/**",
               ":(exclude).RData",
               ":(exclude).Rhistory",
               ":(exclude).Ruserdata")
    if (!biostats_sync_report_git(root, "Stage coursework changes", add)) return(FALSE)

    staged <- git("diff", "--cached", "--quiet")
    if (staged$status == 1L) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
      commit <- git("commit", "-m", paste("Automatic sync:", timestamp))
      if (!biostats_sync_report_git(root, "Commit local changes", commit)) return(FALSE)
      biostats_sync_write(root, "Committed local coursework.")
    } else if (staged$status != 0L) {
      biostats_sync_report_git(root, "Check staged changes", staged)
      return(FALSE)
    }
  }

  for (attempt in seq_len(3L)) {
    fetch <- git("fetch", "--no-tags", "origin", "master")
    if (!biostats_sync_report_git(root, "Fetch origin/master", fetch)) {
      biostats_sync_write(root, "No push attempted because fetch failed.")
      return(FALSE)
    }

    counts <- git("rev-list", "--left-right", "--count", "HEAD...FETCH_HEAD")
    if (!biostats_sync_report_git(root, "Compare local and remote history", counts) ||
        length(counts$output) != 1L) return(FALSE)
    values <- strsplit(trimws(counts$output[[1]]), "[[:space:]]+")[[1]]
    if (length(values) != 2L) {
      biostats_sync_write(root, "ERROR: could not read local/remote commit counts.")
      return(FALSE)
    }
    ahead <- as.integer(values[[1]])
    behind <- as.integer(values[[2]])

    if (behind > 0L) {
      if (ahead > 0L) {
        integrate <- git("rebase", "FETCH_HEAD")
      } else {
        integrate <- git("merge", "--ff-only", "FETCH_HEAD")
      }
      if (!biostats_sync_report_git(root, "Integrate incoming commits", integrate)) {
        aborted <- git("rebase", "--abort")
        if (aborted$status == 0L) {
          biostats_sync_write(root,
            "CONFLICT: incoming edits overlap local work. Rebase was aborted; local commits and files are preserved. Resolve manually, then sync will resume.")
        } else {
          biostats_sync_write(root,
            "INTEGRATION STOPPED: no automatic reset or force-push was attempted. Local work remains for manual recovery.")
        }
        return(FALSE)
      }
      if (behind > 0L) biostats_sync_write(root, "Integrated incoming origin/master commits.")
    }

    ahead_now <- git("rev-list", "--count", "FETCH_HEAD..HEAD")
    if (!biostats_sync_report_git(root, "Check commits to push", ahead_now) ||
        length(ahead_now$output) != 1L) return(FALSE)
    if (as.integer(ahead_now$output[[1]]) == 0L) return(TRUE)

    push <- git("push", "origin", "HEAD:master")
    if (push$status == 0L) {
      biostats_sync_report_git(root, "Push to origin/master", push)
      biostats_sync_write(root, "PUSH SUCCEEDED.")
      return(TRUE)
    }
    biostats_sync_report_git(root, "Push to origin/master", push)
    if (attempt < 3L) {
      biostats_sync_write(root, "Push failed; fetching again before a safe retry.")
    }
  }

  biostats_sync_write(root,
    "PUSH FAILED after 3 attempts. Local commits remain on this computer; check authentication, network access, and the log.")
  FALSE
}

biostats_sync_run_worker <- function(root = biostats_sync_root()) {
  lock <- biostats_sync_lock(root)
  if (!dir.create(lock, showWarnings = FALSE)) {
    biostats_sync_write(root, "Another watcher owns the sync lock; exiting.")
    return(invisible(FALSE))
  }
  on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(as.character(Sys.getpid()), file.path(lock, "pid"))
  Sys.setenv(
    GIT_AUTHOR_NAME = "Joseph Goldenberg",
    GIT_AUTHOR_EMAIL = "ijosephgold@gmail.com",
    GIT_COMMITTER_NAME = "Joseph Goldenberg",
    GIT_COMMITTER_EMAIL = "ijosephgold@gmail.com",
    GIT_TERMINAL_PROMPT = "0",
    GCM_INTERACTIVE = "Never"
  )
  biostats_sync_write(root, "Background watcher started; checking every 60 seconds.")

  repeat {
    if (file.exists(biostats_sync_stop(root))) {
      biostats_sync_write(root, "Stop file detected; watcher stopped.")
      break
    }
    tryCatch(
      biostats_sync_once(root),
      error = function(e) biostats_sync_write(root, "UNEXPECTED ERROR:", conditionMessage(e))
    )
    Sys.sleep(60)
  }
  invisible(TRUE)
}

if ("--worker" %in% commandArgs(trailingOnly = TRUE)) {
  biostats_sync_run_worker()
} else {
  biostats_sync_start()
}
