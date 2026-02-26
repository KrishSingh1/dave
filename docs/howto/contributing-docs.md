# Contributing to the Documentation

This page explains how to add or update pages on this documentation site.

## Overview of the Workflow

1. Write or edit a Markdown file under `docs/`
2. Preview it locally (in VSCode or with `mkdocs serve`)
3. Commit and push to a branch
4. Open a Merge Request on GitLab
5. Once the MR is merged, the site rebuilds automatically

## Writing in Markdown

Documentation pages are plain-text files with a `.md` extension. Markdown is a
lightweight formatting syntax that is readable as-is and renders cleanly in the
browser.

Prefer Markdown over Word for documentation:

- It works with version control — diffs, history, and code review all work
  naturally
- It renders consistently for every reader without font or layout drift
- Equations, code blocks, tables, and cross-links are all first-class features
- No proprietary file format; any text editor works

### Basic Formatting

```markdown
# Heading 1
## Heading 2

**bold**, _italic_, `inline code`

- bullet item
- another item

1. numbered item
2. another item

[link text](other-page.md)
[external link](https://example.com)

![alt text](../path/to/image.png)
```

For mathematics, see [Writing Equations](markdown-math.md).

### Code Blocks

Fenced code blocks use triple backticks and an optional language identifier
for syntax highlighting:

````markdown
```python
def greet(name):
    return f"Hello, {name}"
```
````

### Admonitions

Highlighted callout boxes use the `!!!` syntax:

```markdown
!!! note
    This is a note.

!!! warning
    This is a warning.

!!! tip
    This is a tip.
```

Supported types: `note`, `warning`, `tip`, `info`, `danger`, `success`.


### Tables

```markdown
| Column A | Column B | Column C |
|----------|----------|----------|
| value 1  | value 2  | value 3  |
```

## Previewing Changes

### Option 1: VSCode Markdown Preview

VSCode has a built-in Markdown previewer. Open the `.md` file and press
**Ctrl+Shift+V** (Windows/Linux) or **Cmd+Shift+V** (macOS) to open a
live preview panel beside the editor.

This is the quickest way to check formatting while writing, but it does not
render MkDocs-specific features (admonitions, math, or the nav sidebar).

### Option 2: Local Site Preview with `mkdocs serve`

For a full preview that matches what the deployed site looks like, run the
development server from the repository root:

```bash
pip install mkdocs pymdown-extensions
mkdocs serve --livereload
```

Then open [http://127.0.0.1:8000](http://127.0.0.1:8000) in a browser.
The server watches for file changes and reloads automatically — edit and save
a `.md` file to see the result immediately without restarting.

!!! note
    `mkdocs serve` requires Python. On Windows, use a terminal where Python
    is on `PATH` (e.g. the Anaconda Prompt or a terminal opened from VS Code
    with a Python environment active).

## Alternative: Convert a Word Document

If content is already written in a Word document (`.docx`), the `doc2md`
script in the `docs/` directory can convert it to Markdown using
[pandoc](https://pandoc.org/) (you will first need to install pandoc):

```bash
docs/doc2md input.docx # Creates input.md
```

This handles basic formatting, tables, and equation numbering. Review the
output carefully after conversion — Word's equation objects do not translate
perfectly, and headings, lists, and figure captions often need manual cleanup.

Writing directly in Markdown from the start is generally less work than
converting from Word.

## Adding a New Page

1. Create the `.md` file under `docs/` (e.g. `docs/howto/my-topic.md`).
2. Add it to the `nav:` section of `mkdocs.yml`:
   ```yaml
   - How-to:
       - My Topic: howto/my-topic.md
   ```
3. Verify it appears correctly with `mkdocs serve`.

Pages not listed under `nav:` are still built but will not appear in the
sidebar navigation.

## Committing and Raising a Merge Request

Once you are happy with your changes:

1. Commit the changed files to a branch (see
   [Version Control via SourceTree](vcs-sourcetree.md) or
   [Version Control via CLI](vcs-cli.md))
2. Push the branch to GitLab
3. Open a Merge Request targeting the `lpjguess_dave` branch
4. The CI pipeline runs a compile check; the MR can be merged once checks pass
5. Merging to `lpjguess_dave` triggers the Pages deployment pipeline, which
   rebuilds and publishes the site automatically within a few minutes

No manual deployment step is needed.
