# Common Failures

This page lists common failures, and what to do about them. For the local
VSCode workflow, see [Running and Debugging via VSCode](debugging-vscode.md).

## Start with `guess.log`

When a run fails **or produces an unexpected result**, open `guess.log` first.
For an OzFlux site run, it is in the working directory beside the site
instruction file, for example
`benchmarks/ozflux/CumberlandPlain/guess.log`, **not** in `out/`. Search it for
`FATAL`, `ERROR`, `WARN`, the grid-cell coordinates, the affected stand or PFT,
and the year when the problem appeared. Read the surrounding messages before
changing inputs. The terminal may show a short error, but the log usually has
more context.

Set `verbosity` in the `.ins` file to control how much the log records:

```ins
verbosity 3
```

Use **at least `3`** when investigating a run: this includes informational
messages such as treeinit cohorts and land-cover input decisions. Set `4` for
diagnostic messages or `5` for more detailed debug messages, then rerun and
recheck the log. Higher levels can make large logs. The shared OzFlux
`settings.ins` already sets `verbosity 3`; a site-specific setting can override
it when imported afterward. For parallel runs, each worker writes its own
`runN/guess.log`; see [Running in Parallel](run-parallel.md). On Gadi, see the
[job and worker log guide](gadi-troubleshooting.md).

## First: Did the Run Complete?

When the model finishes normally, the terminal ends with text like:

```text
Simulation completed in 15.05 seconds with 0 errors and 0 warnings.
See guess.log for full details

Finished
```

The time and counts will differ. A completed run can have warnings, and some
model configurations may report non-fatal errors; `Simulation completed` means
the model finished the simulation successfully and wrote outputs.

If the terminal instead contains:

```text
FATAL ERROR:
...
Simulation terminated.
```

the run failed. Note the first line after `FATAL ERROR:`, then open
`guess.log` for the surrounding messages. In VSCode on Windows or Linux, open
the **Terminal** panel at the bottom of the window and scroll up if needed. On
macOS, open **Debug Console** from the panel dropdown instead.

Whether the run completed or failed, check the log in the site's directory.
For Cumberland Plain it is:

```text
benchmarks/ozflux/CumberlandPlain/guess.log
```

Open this file in VSCode, then use `Ctrl+F` (Windows/Linux) or `Cmd+F` (macOS)
to search for `FATAL`, `ERROR`, `WARN`, or a distinctive word from the terminal
message. See [Running and Debugging via VSCode](debugging-vscode.md) for where
outputs and logs are created. Where supported by the build platform, the log also includes a
stack trace after a fatal error; this is mainly useful when reporting the
problem to a developer.

## Specific Common Problems

### `Couldn't find input module nc`

Typical terminal text:

```text
Couldn't find input module nc
```

The executable was built without NetCDF support. The OzFlux simulations use
NetCDF climate forcing, and `nc` is the DAVE input module supplied only when
NetCDF was found during configuration.

Install NetCDF using the setup instructions for your platform, then configure
and build again from the repository root. Do not only rerun the build command:
the configuration step must be repeated so that CMake can find NetCDF.

- [Windows setup](INSTALL.windows.md)
- [macOS setup](INSTALL.macos.md)
- [Linux setup](INSTALL.linux.md)

```bash
cmake -B build -S .
cmake --build build
```

## `NetCDF error: ... Cannot open NetCDF file: ...`

Typical terminal or log text:

```text
NetCDF error: No such file or directory: Cannot open NetCDF file: ../forcing/met/CumberlandPlain.nc
```

The path at the end is the file the model could not find. For the standard
Cumberland Plain run, check that this exact file exists:

```text
benchmarks/ozflux/forcing/met/CumberlandPlain.nc
```

If it is absent, download the full OzFlux forcing directory and place it as
`benchmarks/ozflux/forcing/`; the guide explains the expected location and
transfer options: [HIE-Storage](hie-storage.md).

If the file exists but the error still names a path beginning `../forcing/`,
run the model from the site directory. The command-line steps are in
[Running Serially](run-serial.md).
VSCode's supplied run configurations choose the site directory automatically;
see [Running and Debugging via VSCode](debugging-vscode.md).

## `Error: could not open ... for input`

Typical terminal or log text:

```text
Error: could not open CumberlandPlain.ins for input
```

or:

```text
Failed to open gridlist file: gridlist.txt
```

The model cannot find an instruction or supporting input file. First check the
name in the message for a typo. Then check where the command was run: for an
OzFlux site, the current directory should be that site's folder, for example:

```text
benchmarks/ozflux/CumberlandPlain/
```

The site instruction files deliberately use relative paths. [Running
Serially](run-serial.md) explains this convention, and [OzFlux Site
Simulations](ozflux-simulations.md) describes the files that a site imports.

## `Error: Path does not match any files: ...`

Typical terminal or log text:

```text
Error: Path does not match any files: '../forcing/met/CumberlandPlain.nc'
```

This has the same practical meaning as a missing input file: inspect the path
in the message, confirm that the file is present, and make sure the model was
started from the site directory. Follow the checks in [`NetCDF error: ...
Cannot open NetCDF file: ...`](#netcdf-error-cannot-open-netcdf-file).

## `cmake: command not found`, `cmake is not recognised`, or a Compiler Is Missing

Typical terminal text includes one of:

```text
cmake: command not found
'cmake' is not recognized as an internal or external command
gcc: command not found
```

The build tools are not installed, or the terminal/VSCode was opened before
they were added to the system path. Return to the appropriate platform setup
guide, complete its dependency instructions, then close and reopen the
terminal and VSCode before trying again.

- [Windows setup](INSTALL.windows.md)
- [macOS setup](INSTALL.macos.md)
- [Linux setup](INSTALL.linux.md)

## The Run Configuration Is Missing in VSCode

Symptoms:

- The Run and Debug panel does not offer an OzFlux site such as
  `CumberlandPlain`.
- Pressing `F5` does not build or run an OzFlux simulation.

Open the repository root folder in VSCode, rather than `modules/`,
`benchmarks/`, or an individual site folder. In the Explorer panel, the root
should contain `benchmarks/`, `docs/`, `framework/`, `modules/`, and `.vscode/`.

The tracked `.vscode/launch.json` and `.vscode/tasks.json` files provide the run
configurations. If `.vscode/` is absent, update to the expected branch or ask
Drew to check the checkout. For the normal VSCode workflow, see [Running and
Debugging via VSCode](debugging-vscode.md).

## `Simulation completed ...` but `out/` Is Empty or a File Is Missing

First confirm that you are looking in the site's own output directory, for
example:

```text
benchmarks/ozflux/CumberlandPlain/out/
```

The available outputs are configured in `benchmarks/ozflux/outputs.ins`.
Lines beginning with `!` are disabled. For example:

```text
file_dave_gpp "dave_gpp.out"
!file_dave_gpp "dave_gpp.out"
```

The first line writes `dave_gpp.out`; the second does not. See [OzFlux Site
Simulations](ozflux-simulations.md) for the shared configuration files, and
[Running and Debugging via VSCode](debugging-vscode.md) for the usual output
locations.

## The Run Is Still Printing Progress

If the terminal continues to update a percentage complete, the model is still
running. If `spinup_output` is enabeld, the outputs will include spinup years and may take noticeably longer
in a Debug build. Wait for either the completion text described above or a
fatal error before treating the run as finished.

For a different way to start a local run, see [Running Serially](run-serial.md).
The guides for PBS jobs on Gadi are separate and are not needed for local site
runs.

## SourceTree Shows Many Changed Output Files

Files under `build/`, `benchmarks/ozflux/<SiteName>/out/`, and
`benchmarks/ozflux/<SiteName>/state/` are generated locally. Do not stage them
unless a maintainer has specifically asked for them. The [Repository
Map](repository-map.md) explains which directories contain generated files and
which contain the model inputs or source code.

## After Updating, Files or Run Configurations Are Missing

First check that SourceTree shows the expected branch. If you use the command
line, the following shows the current branch and whether local changes need
attention:

```bash
git status
```

Then fetch and pull the expected branch in SourceTree. The project workflow is
documented in [Version Control via SourceTree](sourcetree.md) and [Version
Control via CLI](git-cli.md). If files are still missing after updating, ask a
maintainer to check whether the checkout is sparse, filtered, or incomplete.
