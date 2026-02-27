# Version Control via the Command Line

This page covers the Git workflow for contributing to LPJ-GUESS DAVE using
the terminal. It assumes familiarity with basic Git concepts (commits,
branches, remotes). For a gentler introduction, see
[Version Control via SourceTree](sourcetree.md).

---

## Clone the Repository

```bash
git clone git@stormbringer4.nateko.lu.se:dave/lpj-guess-dave.git lpj-guess-dave
cd lpj-guess-dave
```

The default branch after cloning is `lpjguess_dave`.

---

## Create a Branch

Always work on a dedicated branch, not directly on `lpjguess_dave`:

```bash
git switch -c my-feature-branch
```

Use a short, descriptive name with hyphens (e.g. `fix-canopy-conductance`,
`add-thome-docs`).

---

## Check Status and Review Changes

```bash
git status          # which files have changed
git diff            # unstaged changes (line-by-line)
git diff --staged   # changes already staged for commit
```

---

## Stage and Commit

Stage specific files you want to include in the next commit:

```bash
git add path/to/changed/file.cpp path/to/file2.h
git add docs/howto/my-page.md
```

Or stage everything in the working tree (use with care):

```bash
git add -A   # stage all changes in the repository, including deletions
git add .    # stage all changes in the current directory and subdirectories
git add -p   # interactive patch staging — recommended over git add .
```

Then commit with a message that describes what changed and why:

```bash
git commit -m "Fix negative VPD crash when RH input is 1.0"
```

Good commit messages use the imperative present tense and explain why the
change was made, not just which files were touched.

---

## Push the Branch

```bash
git push -u origin my-feature-branch
```

The `-u` flag sets the upstream tracking reference so that subsequent
`git push` and `git pull` calls on this branch need no arguments.

---

## Keep the Branch Up to Date

If `lpjguess_dave` has moved forward while you were working, rebase your
branch onto it to avoid conflicts at merge time:

```bash
git fetch origin
git rebase origin/lpjguess_dave
```

If there are conflicts, Git will pause and ask you to resolve them file by
file. After resolving each conflict:

```bash
git add <resolved-file>
git rebase --continue
```

### Alternative: Merge Instead of Rebase

If you prefer not to rewrite your branch history, merge `lpjguess_dave` into
your branch instead:

```bash
git fetch origin
git merge origin/lpjguess_dave
```

If there are conflicts, Git will pause and ask you to resolve them file by
file. After resolving each conflict:

```bash
git add <resolved-file>
git merge --continue
```

If your Git version does not support `git merge --continue`, complete the
merge with a regular commit:

```bash
git commit
```

### Rebase vs Merge (What’s the Difference?)

- `rebase` replays your branch commits on top of the latest
    `origin/lpjguess_dave`, producing a linear history.
- `merge` creates a merge commit that ties both histories together, preserving
    the exact original branch timeline.

Practical guidance:

- Prefer `rebase` for a clean, linear feature branch history before opening an
    MR.
- Prefer `merge` when you want to avoid history rewriting, especially on
    shared branches.
- Avoid rebasing commits that others have already based work on.

---

## Open a Merge Request

After pushing, open the repository in your browser on GitLab. GitLab will
show a banner prompting you to create a Merge Request for the branch you
just pushed. Set **Target branch** to `lpjguess_dave`.

Alternatively, GitLab prints a direct MR creation link in the output of
`git push`:

```
remote: To create a merge request for my-feature-branch, visit:
remote:   https://gitlab.com/.../merge_requests/new?merge_request[source_branch]=my-feature-branch
```

The CI pipeline runs automatically. The MR can be merged once the pipeline
passes and a reviewer has approved it.

!!! note
    Subsequent pushes to the same branch update the open MR automatically —
    no need to open a new one.

---

## Pull Updates from the Remote

To bring in commits merged to `lpjguess_dave` by others:

```bash
git switch lpjguess_dave
git pull
```

---

## Quick Reference

| Task | Command |
|------|---------|
| Clone | `git clone <url>` |
| New branch | `git switch -c <name>` |
| Switch branch | `git switch <name>` |
| Stage file | `git add <file>` |
| Interactive staging | `git add -p` |
| Commit | `git commit -m "<message>"` |
| Push (first time) | `git push -u origin <branch>` |
| Push (subsequent) | `git push` |
| Fetch + rebase | `git fetch origin && git rebase origin/lpjguess_dave` |
| Fetch + merge | `git fetch origin && git merge origin/lpjguess_dave` |
| Pull main branch | `git switch lpjguess_dave && git pull` |
| Check status | `git status` |
| View diff | `git diff` / `git diff --staged` |
| View log | `git log --oneline` |
