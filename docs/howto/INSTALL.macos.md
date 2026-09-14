# macOS Setup Guide

This guide covers how to set up and compile LPJ-GUESS DAVE on macOS using the
Clang/Apple toolchain and VSCode. After completing this setup, follow
[Running and Debugging via VSCode](debugging-vscode.md) to run a site
simulation.

## Software Installation

Many setup steps use Terminal. To open Terminal, press `Cmd+Space`, type
`Terminal`, and press Enter. When this guide says to run a command, paste it
into Terminal and press Enter.

### Xcode Command Line Tools

The Xcode Command Line Tools provide the compiler (clang, which masquerades as
gcc) and other essential build utilities. To check whether they are already
installed, run:

```bash
gcc --version
```

If they are not installed, you will see a dialog prompting you to install them.
You can also trigger this manually:

```bash
xcode-select --install
```

### Homebrew

Homebrew is a package manager for macOS used to install the remaining
dependencies. To check if it is already installed:

```bash
brew --version
```

If it is not installed, follow the instructions at
[brew.sh](https://brew.sh).

### Dependencies

Install CMake and NetCDF via Homebrew:

```bash
brew install cmake netcdf
```

### Visual Studio Code

Download and install VSCode from the
[official website](https://code.visualstudio.com/Download).

Once installed, open the Extensions panel (`Cmd+Shift+X`) and install the
following extensions:

- **C/C++** — language support and IntelliSense
- **CodeLLDB** — debugger integration (uses lldb under the hood)
- **Command Variable** — used by the run/debug launch configurations to prompt
  for input

## Obtaining Source Code

Clone the DAVE repository:

```bash
git clone git@stormbringer4.nateko.lu.se:dave/lpj-guess-dave.git
```

Or via HTTPS if you have not set up an SSH key:

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
