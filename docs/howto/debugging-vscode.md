# Running and Debugging via VSCode

The repository includes pre-configured VSCode launch and build tasks in
`.vscode/launch.json` and `.vscode/tasks.json`. These are committed to the
repository, so debugging should work without any additional configuration after
completing the platform setup.

## Before You Start

- VSCode installed with the required extensions (see your platform setup guide:
  [Windows](INSTALL.windows.md) | [macOS](INSTALL.macos.md) |
  [Linux](INSTALL.linux.md))
- Climate forcing data downloaded to `benchmarks/ozflux/forcing/`
  (see [HIE-Storage](hie-storage.md))

For the example below, check that this file exists:

```text
benchmarks/ozflux/forcing/met/CumberlandPlain.nc
```

If it does not, download the full OzFlux forcing directory before continuing.
The [OzFlux Site Simulations](ozflux-simulations.md) page explains the input
file structure.

## Run Your First Site

1. Open the repository root folder in VSCode (`File → Open Folder`). Always
   open the root directory — not a subdirectory such as `modules/` or
   `benchmarks/`.

2. Open the Run and Debug panel (`Ctrl+Shift+D` on Windows/Linux,
   `Cmd+Shift+D` on macOS), or click the debug icon in the left sidebar, or
   click View -> Run.

3. Select `CumberlandPlain` from the run-configuration dropdown at the top of
   the panel. You can select a different OzFlux site later.

4. Press `F5` (or the green play button) to start.

The pre-launch build task will recompile the model if any source files have
changed since the last build. Compilation output appears in the Terminal panel.
On Windows and Linux, the model then runs in that terminal, with progress
reported there. On macOS, the debugged model's standard output and standard
error appear in the **Debug Console** panel instead; select **Debug Console**
from the panel dropdown at the bottom of VSCode.

When the run reaches the end, it prints text like:

```text
Simulation completed in 15.05 seconds with 0 errors and 0 warnings.
See guess.log for full details
```

The time and counts will differ. `FATAL ERROR:` followed by `Simulation
terminated.` means the run failed; see [Common Failures](common-failures.md).

## Check the Results

For Cumberland Plain, outputs are written to:

```text
benchmarks/ozflux/CumberlandPlain/out/
```

The full log is:

```text
benchmarks/ozflux/CumberlandPlain/guess.log
```

Open the log first if the run failed, reported warnings or errors, **or simply
gave an unexpected result**. Search for the relevant year, PFT, stand, and
warning or error. `verbosity` in the `.ins` file controls the detail: use at
least `3`, or `4`/`5` for more detail. See
[Start with `guess.log`](common-failures.md#start-with-guesslog). The exact output files depend
on the shared `benchmarks/ozflux/outputs.ins` configuration; see [OzFlux Site
Simulations](ozflux-simulations.md) before changing it. Do not commit files in
`out/` unless a maintainer has asked for them.

To run another site, select its configuration from the same dropdown and press
`F5`. For command-line runs, see [Running Serially](run-serial.md).

## Debugger Notes

- **Linux** uses `gdb`; **macOS** uses `lldb`; **Windows** uses `gdb` via
  MSYS2. The launch configuration selects the correct debugger automatically
  based on platform.
- To set a breakpoint: click in the left margin next to a line number in the
  editor, then press `F5`. Execution will pause at that line.
- The debug console at the bottom of the screen accepts GDB/LLDB commands
  directly.
