# Running Serially (CLI)

This page explains how to run a simulation from the command line without an IDE
or a debugger.

## Basic Syntax

```bash
/path/to/build/guess -input <module> <instruction-file>
```

For DAVE simulations, the input module is always `nc` (see
[Climate Input](climate-input.md)):

```bash
cd /path/to/repository
cd benchmarks/ozflux/CumberlandPlain
../../../build/guess -input nc CumberlandPlain.ins
```

In general, the working directory is recommended to be the same as the
instruction file. This is because instruction files can contain relative paths,
which are resolved relative to the working directory, but when writing
instruction files, one tends to think in terms of the instruction file itself.
Therefore, it's easiest to just run from the same directory as the instruction
file.

## OzFlux Sites

Each OzFlux site has its own instruction file at
`benchmarks/ozflux/<SiteName>/<SiteName>.ins`. See
[OzFlux Site Simulations](ozflux-simulations.md) for details of the input
file structure.

## Useful Options

| Flag | Description |
|---|---|
| `-input <module>` | Input module to use (typically `nc` when running DAVE, often `cf` when running upstream LPJ-GUESS with daily input) |
| `-parallel` | Enable MPI parallelism (requires MPI build; see [parallel runs](run-parallel.md)) |

## Compiled Binary Location

After a standard CMake build the binary is at `build/guess` (Unix) or
`build\guesscmd.exe` (Windows).
