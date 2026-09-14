# Linux Setup Guide

This guide covers how to set up and compile LPJ-GUESS DAVE on Linux. After
completing this setup, follow [Running and Debugging via VSCode](debugging-vscode.md)
to run a site simulation.

## Dependencies

On Debian/Ubuntu:

```bash
sudo apt install cmake g++ libnetcdf-dev
```

On other distributions, install the equivalent packages via your native package
manager. You need:

- A C++ compiler (GCC or Clang)
- CMake
- NetCDF-C development headers and library

## Visual Studio Code (optional)

VSCode is recommended for interactive debugging (though not strictly required).
Download it from the
[official website](https://code.visualstudio.com/Download), then install the
following extensions (`Ctrl+Shift+X`):

- **C/C++** — language support and IntelliSense
- **CodeLLDB** — debugger integration
- **Command Variable** — used by the run/debug launch configurations to prompt
  for input

## Obtaining Source Code

Clone the DAVE repository:

```bash
git clone git@stormbringer4.nateko.lu.se:dave/lpj-guess-dave.git
```

Or via HTTPS:

```bash
git clone https://stormbringer4.nateko.lu.se/dave/lpj-guess-dave.git
```

## Compilation

From the repository root:

```bash
cmake -B build -S .
cmake --build build
```

The compiled binary will be at `build/guess`.

## Next Step: Run the Model

Download the OzFlux forcing data from [HIE-Storage](hie-storage.md), then follow
[Running and Debugging via VSCode](debugging-vscode.md) to run the
`CumberlandPlain` site and check the output files.
