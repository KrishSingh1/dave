# Debugging via VSCode

The repository includes pre-configured VSCode launch and build tasks in
`.vscode/launch.json` and `.vscode/tasks.json`. These are committed to the
repository, so debugging should work without any additional configuration after
completing the platform setup.

## Prerequisites

- VSCode installed with the required extensions (see your platform setup guide:
  [Windows](INSTALL.windows.md) | [macOS](INSTALL.macos.md) |
  [Linux](INSTALL.linux.md))
- Climate forcing data downloaded to `benchmarks/ozflux/forcing/`
  (see [HIE-Storage](hie-storage.md))

## Running a Site

1. Open the repository root folder in VSCode (`File → Open Folder`). Always
   open the root directory — not a subdirectory such as `modules/` or
   `benchmarks/`.

2. Open the Run and Debug panel (`Ctrl+Shift+D` on Windows/Linux,
   `Cmd+Shift+D` on macOS), or click the debug icon in the left sidebar, or
   click View -> Run.

3. Select a site from the **Run Site** dropdown at the top of the panel.

4. Press `F5` (or the green play button) to start.

The pre-launch build task will recompile the model if any source files have
changed since the last build. Compilation output appears in the Terminal panel.
The model then runs in the same terminal, with progress reported to stdout.

## Outputs

Outputs are written to `benchmarks/ozflux/<SiteName>/out/` by default.

## Debugger Notes

- **Linux** uses `gdb`; **macOS** uses `lldb`; **Windows** uses `gdb` via
  MSYS2. The launch configuration selects the correct debugger automatically
  based on platform.
- To set a breakpoint: click in the left margin next to a line number in the
  editor, then press `F5`. Execution will pause at that line.
- The debug console at the bottom of the screen accepts GDB/LLDB commands
  directly.
