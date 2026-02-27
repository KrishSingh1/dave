# Version Control via SourceTree

[SourceTree](https://www.sourcetreeapp.com/) is a free graphical Git client
for Windows and macOS. This page walks through the full workflow for
contributing to LPJ-GUESS DAVE — from cloning the repository to opening a
Merge Request on GitLab.

---

## Key Concepts

Git manages a complete history of changes to the codebase. Understanding a
few terms will make the rest of this guide easier to follow:

**Repository (repo)**
: The project folder plus the complete record of every change ever made to
  it. The authoritative copy lives on GitLab; everyone who works on the code
  has their own local copy.

**Commit**
: A saved snapshot of the repository at a point in time, with a short
  description of what changed and why. Commits are permanent and can always
  be recovered.

**Branch**
: An isolated line of development. The main branch (`lpjguess_dave`) always
  contains tested, reviewed code. You do your work on a separate branch so
  that the main branch is never broken by work in progress.

**Staging area**
: Before a commit is made, you choose exactly which changed files (and even
  which lines within a file) to include. This is called *staging*. It lets
  you make several unrelated edits and then commit them separately with
  meaningful descriptions.

**Remote**
: The shared copy of the repository on GitLab. "Pushing" sends your local
  commits to the remote; "pulling" brings in commits others have pushed.

**Merge Request (MR)**
: A GitLab feature. When your branch is ready, you open an MR to ask for it
  to be reviewed and merged into `lpjguess_dave`. The CI pipeline runs
  automatically and reviewers can leave comments on your changes.

---

## Installing SourceTree

Download and install SourceTree from
[sourcetreeapp.com](https://www.sourcetreeapp.com/). During setup it will
ask you to sign in with an Atlassian account (free).

When prompted to set up a hosting account, choose **GitLab** and
authenticate with your GitLab credentials. This allows SourceTree to push
and pull without asking for your password on every operation.

---

## Cloning the Repository

Cloning creates a local copy of the repository on your machine.

1. In SourceTree, open **File → New...** and select the **Clone** tab.
2. Paste the repository URL into the **Source Path** field:
   ```
   <repository-url>
   ```
3. Set **Destination Path** to a folder on your machine
   (e.g. `C:\Users\you\code\lpj-guess-dave`).
4. Leave **Name** as the default, then click **Clone**.

![SourceTree clone dialog](img/sourcetree/clone.png)

SourceTree will download the full repository. This may take a minute or two
on the first clone.

---

## Overview of the SourceTree Window

Once cloned, SourceTree opens the repository view.

![SourceTree main window](img/sourcetree/main.png)

The key areas are:

- **Toolbar** — buttons for the most common actions: Commit, Pull, Push,
  Fetch, Branch.
- **Left sidebar** — lists your local branches, remote branches, and tags.
  The branch you currently have checked out is shown in bold.
- **Commit graph** (centre) — shows the history of commits as a timeline.
  Each dot is a commit; branches and merges are visible as lines.
- **Working Copy** entry (top of the history list) — click this to see
  your uncommitted changes.

---

## Creating a Branch

Always do new work on its own branch, never directly on `lpjguess_dave`.

1. Click the **Branch** button in the toolbar.
2. Enter a short, descriptive name in the **New Branch** field
   (e.g. `fix-canopy-conductance` or `add-thome-docs`). Use hyphens, not
   spaces.
3. Make sure **Checkout New Branch** is ticked.
4. Click **Create Branch**.

![SourceTree new branch dialog](img/sourcetree/new-branch.png)

The sidebar will now show your new branch in bold, meaning it is the active
branch. Any commits you make will go onto this branch.

---

## Making and Staging Changes

Edit files in the repository folder using whatever tools you normally use
(an IDE, VSCode, a text editor). When you return to SourceTree you will see
your changes listed.

Click **Working Copy** at the top of the commit graph to open the staging
view:

![SourceTree staging area](img/sourcetree/working-copy.png)

The panel has two sections:

- **Unstaged files** (top) — files you have changed but not yet selected
  for the next commit.
- **Staged files** (bottom) — files that will be included in the next
  commit.

To stage a file, click the checkbox to the left of its name in the
**Unstaged files** panel. It moves to **Staged files**. You can stage all
changed files at once with the **Stage All** button.

The diff viewer on the right shows exactly what changed in the selected
file — lines added in green, lines removed in red.

!!! tip
    Stage only the files related to one logical change. If you made two
    unrelated edits, stage and commit them separately with different
    descriptions. This makes the history easier to read and simplifies
    reviewing.

---

## Committing

Once you have staged the files you want, write a commit message in the text
box at the bottom of the staging view.

A good commit message:

- Is a single sentence in the imperative present tense:
  *"Fix negative VPD crash when RH input is 1.0"*
  *"Add temperature acclimation documentation"*
- Describes **why** the change was made, not just what files were touched.

Click **Commit** to save the snapshot. The commit appears immediately in the
graph — it is saved locally, but not yet on GitLab.

---

## Pushing to GitLab

"Pushing" sends your local commits to GitLab so others can see them and so
that the CI pipeline can run.

1. Click the **Push** button in the toolbar.
2. In the dialog, make sure your branch is ticked.
3. Click **Push**.

![SourceTree push dialog](img/sourcetree/push.png)

If this is the first time you have pushed this branch, GitLab will create
a remote copy of it automatically.

---

## Opening a Merge Request

Once your branch is pushed:

1. Open the repository on GitLab in your browser.
2. GitLab will usually show a banner: **"You pushed a branch — create a
   Merge Request."** Click the **Create Merge Request** button.
3. Fill in the title and description. The description should explain what
   the change does and how to test it.
4. Set **Target branch** to `lpjguess_dave`.
5. Submit the MR.

The CI pipeline will run automatically. Your MR can be merged once the
pipeline passes and a reviewer has approved it.

!!! note
    You do not need to open a new MR every time you push. If you push
    additional commits to the same branch, they appear in the already-open
    MR automatically.
