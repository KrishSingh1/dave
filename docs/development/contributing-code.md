# Contributing Code (Policy)

This page defines governance and quality expectations for code contributions.
For step-by-step contributor workflow, see
[Contributing Code (Workflow)](../howto/contributing-code.md).

For command-level Git instructions, see:

- [Version Control via CLI](../howto/git-cli.md)
- [Version Control via SourceTree](../howto/sourcetree.md)

---

## 1. Roles and Permissions

### Maintainers (ie Drew)

- Maintain this fork and its protected branches.
- Sync changes from the upstream LPJ-GUESS repository.
- Review and merge Merge Requests (MRs).
- Curate the overall structure of long-lived development branches.

### Developers (ie postdocs)

- Implement new features and fixes.
- Work primarily on long-lived personal branches in this team repository.
- Open MRs from long-lived branches back to `lpjguess_dave`.

### Reporters (ie students)

- Run the model and test results.
- Report bugs and issues.
- Review documentation and suggest improvements.
- Optional: contribute code via branches (with guidance from a maintainer).

---

## 2. Branch and Repository Policy

### Protected branches

- `lpjguess_dave`: primary, stable branch for team work. All merges require MR and review.
- `trunk`: tracking branch for upstream LPJ-GUESS. Maintainers only.

### Long-lived personal branches

Naming convention:

- `lpjguess_dave_<name>`
- `lpjguess_dave_<name>_<short-topic>` (optional qualifier)

Examples:

- `lpjguess_dave_clare`
- `lpjguess_dave_juergen_adaptive`

Policy:

- Long-lived branches live in this team repository (not personal forks).
- They remain the primary working branch for each developer.
- They are not deleted after merges to `lpjguess_dave`.

### Short-lived feature branches

- Use for risky or isolated work.
- Merge into either the relevant long-lived branch or directly into `lpjguess_dave`.
- Delete after merge.

---

## 3. Upstream Sync Policy (Maintainers)

Maintainers periodically sync from upstream LPJ-GUESS to `trunk`, then open
an MR from `trunk` to `lpjguess_dave`.

After this integration, developers should update their own long-lived branches.

---

## 4. MR Quality Expectations

All MRs targeting `lpjguess_dave` should include:

- What changed.
- Why it changed (scientific or technical motivation).
- How it was tested (runs/configs/checks).
- Expected impact on model outputs, where relevant.

Keep commits and MRs focused; prefer smaller reviewable changes.

---

## 5. Issues and Bug Reporting

Use GitLab Issues for:

- Bugs
- Performance problems
- Documentation gaps
- Feature requests

Include:

- Clear description
- Reproduction steps
- Expected vs actual behaviour
- Relevant input files/logs
- Branch context (`lpjguess_dave`, `lpjguess_dave_<name>`, etc.)

---

## 6. Code and Scientific Documentation Practices

- Write clear, readable code.
- Comment non-obvious logic or scientific assumptions.
- Document important equations and references in code comments or Markdown.
- Highlight output-impacting changes in MR descriptions.

Before merge to `lpjguess_dave`, ensure code compiles/runs and outputs pass
sanity checks appropriate to the change.

