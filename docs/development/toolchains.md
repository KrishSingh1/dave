# Dependencies and Toolchains

This page describes the build system, compilers, and library dependencies used
to compile LPJ-GUESS DAVE. It is intended as an orientation for developers
rather than a step-by-step setup reference — see the
[platform setup guides](#platform-setup-guides) for the latter.

## Build System

DAVE (and LPJ-GUESS more broadly) uses **CMake**. CMake generates native build
files (Makefiles on Unix-like systems, Ninja files, or IDE project files) from a
platform-independent description. All compiler and dependency detection happens
at the CMake configure step.

All external library dependencies are detected automatically and are
**optional at the CMake level** — the project will configure and compile
without any of them, but certain features will be disabled. The CMake output
will report which dependencies were found.

## Compilers

The following compilers are supported:

| Platform | Compiler | Notes |
|---|---|---|
| Linux | GCC | Primary development target |
| macOS | Clang (via Xcode) | Exposed as `gcc` by the Xcode Command Line Tools |
| Windows | MinGW-w64 (GCC) | Installed via MSYS2; the recommended Windows toolchain |
| Windows | MSVC | Technically supported but not extensively tested |

## Library Dependencies

### NetCDF (practically required)

The NetCDF-C library is required to read NetCDF input files. Since all standard
DAVE simulation workflows (OzFlux site-level and beyond) use NetCDF-format
climate forcing data, the model is not useful in practice without this
dependency.

### ZLIB (optional)

If ZLIB is found, output file compression is enabled automatically. This is
a minor convenience feature and has no effect on science outputs.

### MPI (optional)

If an MPI implementation is found on Unix-like systems, MPI support is compiled
in, enabling parallel job submission workflows. MPI is not needed for
single-core development or debugging runs.

## Platform Setup Guides

For installation instructions on your specific platform:

- [Windows](../howto/INSTALL.windows.md)
- [macOS](../howto/INSTALL.macos.md)
- [Linux](../howto/INSTALL.linux.md)
