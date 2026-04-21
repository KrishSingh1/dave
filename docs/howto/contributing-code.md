# Contributing Code (Workflow)

This page covers the practical day-to-day workflow for contributing code.
For branch policy, roles, and maintainer responsibilities, see
[Contributing Code (Policy)](../development/contributing-code.md).

If you need UI/terminal-specific Git guidance, see:

- [Version Control via SourceTree](sourcetree.md)
- [Version Control via CLI](git-cli.md)

---

## 1. Getting Started (First-Time Setup)

Clone the DAVE team repository and confirm you can authenticate and push.
Use either the CLI or SourceTree guides above for exact commands.

If you previously worked on the **old/core LPJ-GUESS repo**, see
[Migrating from the Old Upstream LPJ-GUESS Repository](migrate-from-upstream-repo.md).

---

## 2. Primary Workflow: Long-Lived Personal Branches

This is the **main workflow** for most developers (postdocs).

### 2.1. Create or update your long-lived branch

Make sure your local `lpjguess_dave` is current, then create or switch to your
long-lived branch (`lpjguess_dave_<your_name>`).

Periodically bring changes from `lpjguess_dave` into your long-lived branch.
If you are unsure when to merge vs rebase, ask a maintainer.

---

### 2.2. Make and commit changes

Work as usual in your `lpjguess_dave_<your_name>` branch and commit focused,
logical changes.

Keep commits reasonably small and focused where possible.

Push regularly so your branch is backed up and ready for review.

---

### 2.3. Creating a Merge Request from your dev branch to `lpjguess_dave`

Periodically, when your branch is in a good state to integrate:

1. Ensure your branch is up to date.
2. In GitLab, open a **Merge Request**:
   - Go to the gitlab repository -> code -> merge requests -> new (or use [this link](https://stormbringer4.nateko.lu.se/dave/lpj-guess-dave/-/merge_requests/new))
   - Source branch: `dave/lpj-guess-dave` : `lpjguess_dave_<your_name>`
   - Target branch: `dave/lpj-guess-dave` : `lpjguess_dave`

   NOTE: The repository (first dropdown) in the target branch defaults to
   `lpj-guess-developers/LPJ-GUESS`. This needs to be changed to
   `dave/lpj-guess-dave`; otherwise, the MR will be submitted against the
   upstream LPJ-GUESS repository.
3. Fill out the MR description:
   - What the change does
   - Why it’s needed (scientific or technical motivation)
   - How you tested it (runs / configs)
   - Any expected impacts on model outputs

When the MR is approved and tests pass, a maintainer will merge it.

Your `lpjguess_dave_<your_name>` branch remains your main working branch;
you **do not delete** it after merges.

---

## 3. Secondary Workflow: Short-Lived Feature or Experiment Branches

(For advanced or experimental work.)

### 3.1. Create a short-lived branch

Create a short-lived `feature/<topic>` branch from either your long-lived
branch or `lpjguess_dave`.

Work and commit normally, then open an MR:

- **feature → lpjguess_dave_<your_name>**, or
- **feature → lpjguess_dave** (if ready for direct integration)

After merge, delete the feature branch (GitLab can automate this).

---

## 4. Where to find Git commands

This page intentionally focuses on workflow and review gates.
For concrete command sequences, use:

- [Version Control via CLI](git-cli.md)
- [Version Control via SourceTree](sourcetree.md)
