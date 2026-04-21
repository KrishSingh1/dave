# Migrating from the Old Upstream LPJ-GUESS Repository

This is for people who previously worked on the **original/core LPJ-GUESS repo**
and want to repoint their existing local clone to this new Australian team fork.

If you’re cloning fresh, skip this section.

---

## 1. Check current remote

```bash
git remote -v
```

You’ll likely see something like:

```text
origin git@stormbringer4.nateko.lu.se:lpj-guess-developers/LPJ-GUESS.git (fetch)
origin git@stormbringer4.nateko.lu.se:lpj-guess-developers/LPJ-GUESS.git (push)
```

---

## 2. Commit any local work

Make sure no work is lost:

```bash
git status
```

If you see uncommitted changes:

```bash
git add <files>
git commit -m "Save work before migration"
```

---

## 3. Point origin to the new team fork

Replace the URL below with the actual team fork URL:

```bash
git remote set-url origin git@stormbringer4.nateko.lu.se:dave/lpj-guess-dave.git
```

Verify:

```bash
git remote -v
```

We intentionally do **not** keep a remote pointing to the old upstream, to avoid
accidental pushes.

---

## 4. Push your branches to the new fork

List local branches:

```bash
git branch
```

Push each branch you care about:

```bash
git checkout <branch-name>
git push -u origin <branch-name>
```

---

## 5. Rename your long-lived branch (optional but recommended)

If your main working branch was not named `lpjguess_dave_<your_name>`:

```bash
git checkout <old-branch>
git branch -m lpjguess_dave_<your_name>
git push -u origin lpjguess_dave_<your_name>
git push origin --delete <old-branch>   # optional
```

Ask a maintainer before deleting branches on the remote.

---

## 6. Final checks

- `git fetch`, `git pull`, and `git push` should work normally.
- `git branch -a` should show:
  - `lpjguess_dave`
  - `trunk`
  - your `lpjguess_dave_<name>` branch

At this point, you’re fully migrated and can follow the primary workflow in
[Contributing Code (Workflow)](contributing-code.md).

---

If you have improvements to suggest, please open an Issue titled:

**“Update CONTRIBUTING / migration docs”**

