# Safe automatic synchronization on both Windows laptops

The project startup file launches `run_sync.R` in a separate R process, so the RStudio console remains available. The watcher checks every 60 seconds. It commits coursework locally, fetches `origin/master`, integrates incoming commits, and uses a normal fast-forward push. It never force-pushes or resets local commits.

## Before updating either laptop

1. Save open scripts and data files in RStudio.
2. In File Explorer, make a backup copy of these items from the Biostatistics project folder to a location outside the repository, such as Documents\Biostatistics-session-backup:
   - `.Rproj.user`
   - `.RData`
   - `.Rhistory`
   - `.Ruserdata`
   
   The tracked `.Rproj.user` folder contains RStudio job logs and saved source buffers. Keep the backup until both laptops have been updated and you have confirmed your work is present.
3. Update the repository files on each laptop using GitHub Desktop: select the Biostatistics repository, choose **Fetch origin**, then **Pull origin**. Open the project by double-clicking its `.Rproj` file. The watcher starts automatically.
4. If the watcher reports an incoming-change conflict, do not reset or discard anything. Read `.git/automatic-sync.log`, save your open work, and resolve the conflict in RStudio or GitHub Desktop. The watcher aborts a conflicting rebase and retains the local commits.

## Stop tracking RStudio session state

The ignore rules do not remove files already tracked by Git. After making the backups above on both laptops, open Git Bash in the Biostatistics folder and run:

```sh
git rm -r --cached --ignore-unmatch .Rproj.user .RData .Rhistory .Ruserdata
```

This removes those paths from Git's index while leaving their local copies in place. The automatic watcher will commit the index cleanup and push it. On the other laptop, the incoming deletion may be blocked by changed local session files; make its backup, run the same command there, then let the watcher sync. Do not delete the backups until both laptops show a clean working tree and the cleanup commit is on GitHub.

## Each laptop

- Keep GitHub Desktop installed and signed in to the same GitHub account. Git push uses the Git credential configured for Git on that laptop; the watcher disables terminal prompts so missing/expired credentials are reported in the log instead of hanging.
- Open the same cloned repository through its `.Rproj` file. Do not create a second independent copy or edit the same file on both laptops at the same time.
- To pause the watcher, create an empty file at `.git/automatic-sync.stop`; to resume, delete that file and reopen the project. The log is at `.git/automatic-sync.log`.
- If a push fails, the local commit stays on that laptop. Check the log for authentication, network, or conflict messages, reconnect GitHub Desktop if needed, and leave the local branch intact. Never force-push.
- If you already have multiple local commits or an unfinished Git operation, resolve that before relying on automatic synchronization.

## Checks

The watcher confirms Git command exit codes before reporting a push success. It fetches before each push attempt, retries a changed remote up to three times using fetch and rebase, and aborts a conflicting rebase to preserve local commits. Authentication and network failures are recorded in the per-computer log. RStudio session files are excluded from automatic staging, including when they were previously tracked.
