# LPJ-Guess DAVE

This is the DAVE build of LPJ-Guess.

This version of lpj-guess is configured for debugging with visual studio code,
via gdb or lldb. The following instructions describe the installation of the
required software dependencies and running a simulation of an ozflux site with a
debugger attached.

This branch (dave-baseline) behaves the same as the trunk version, with the
addition of patch-level outputs and some extra code for running the ozflux site-level simulations.

## Dependencies

The following software packages are required in order to compile and debug the
model. Installation instructions for each dependency are given below for MacOS
and Windows. Linux users should just use their native package manager. Users of
other OSs are on their own.

- gcc (compilers)
- Homebrew (a package manager, used to install some of the other tools) (**only required on MacOS**)
- NetCDF libraries
- CMake (build system which calls the compilers)
- Subversion (to get the code)
- Visual Studio Code

### GCC

**Windows**

1. Install [msys2](https://www.msys2.org/)
2. Launch msys2 mingw64 prompt from the start menu
3. Run `pacman -Sy mingw-w64-x86_64-gcc

**MacOS**

Need to first install the xcode dev tools. You can tell if you've done this by attempting to run any gcc command (e.g. `gcc --version`). If it's not installed, you will get output looking like:

```
xcode-select: no developer tools were found at '/Applications/Xcode.ap', requesting install. Choose an option in the dialog to download the command line developer tools.
```

From here, either follow the dialog, or else start the dialog manually with `xcode-select --install`

### Homebrew (MacOS-only)

To check if homebrew is already installed, try running a command such as
`brew --version`

If homebrew is not installed, it may be installed by running something like:

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### Subversion

This is optional (but still useful) if another svn client (e.g. SmartSVN) is
already installed.

**Windows**: Download and run the windows installer from
[their website](https://subversion.apache.org/packages.html#windows).

**MacOS**: `brew install subversion`

### CMake

**Windows**:

Launch msys2 mingw64 prompt from the start menu and run:

`pacman-Sy mingw-w64-x86_64-cmake`

**MacOS**: `brew install cmake`

### NetCDF

Technically this is optional, but it's necessary in order to run simulations
with netcdf input files (such as the ozflux site-level simulations).

**Windows**: Download and run the NetCDF4 installer from
[here](https://docs.unidata.ucar.edu/netcdf-c/current/winbin.html).

**MacOS**: `brew install netcdf`

### Visual Studio Code

This is highly recommended for development and debugging, though ultimately the
model may be compiled and run without this.

To install it, download and run the installer from
[here](https://code.visualstudio.com/Download).

Once vscode is installed, the following extensions must be installed. To install 
an extension, open the extensions tab (Click "View -> Extensions"), search for
the extension using the names provided below, and click "install".

- C/C++
- CodeLLDB
- Command Variable
- SVN (optional, but useful)

## Obtaining the Source Code

The URL of the dave baseline branch is:

`svn://stormbringer.nateko.lu.se/svn/LPJ-GUESS/branches/lpjguess_dave_baseline`

Instructions are provided below for obtaining the source code via SmartSVN
(a graphical svn client) and the svn command-line interface. Other svn clients
such as TortoiseSVN will also work, but usage of these tools is left as an
exercise to the enthusiastic reader.

### Obtaining the source code via SmartSVN

0. Open SmartSVN. The welcome dialog should popup, but if it doesn't, click
   "Project -> Check Out" and proceed to step 2
1. Select "check out project from repository", and click OK
2. Copy the branch URL above into the "Repository" text input, and click Next
3. Ensure that the lpjguess_dave_baseline branch is selected and click Next
4. Into the "Local Directory" input, enter the desired path to the source code
   on your computer. Ensure "Checkout Depth" is set to "Fully Recursive", and
   click Next
5. The default options in the "Project" screen are fine, so click Next unless
   you know what you're doing

### Obtaining the source code via CLI

To obtain the source code via the command-line, run the following command:

```bash
url=svn://stormbringer.nateko.lu.se/svn/LPJ-GUESS/branches/lpjguess_dave_baseline
svn checkout $url /path/to/lpj/guess
```

This will do a checkout into `/path/to/lpj/guess`, so this should be changed to
whatever path is desirable.

## Compilation

### Compilation in Visual Studio Code

When using the model with visual studio code (henceforth vscode), always open
the "main" lpj-guess directory via the "File -> Open Folder" button. The "main"
directory is the one which contains directories such as modules, framework, etc.
Do not open one of the subdirectories (such as modules) - if you do this,
compilation and debugging will not work.

The model may be compiled "manually" by using the Ctrl + Shift + B (Cmd + Shift + B on MacOS) keyboard shortcut. Note however that this is rarely necessary,
because the model will be automatically recompiled whenever it is launched in
vscode with the debugger (see the debugging section below).

### Compilation from the CLI

To compile the model from the command line, run these commands, susbsituting the
appropriate path to the lpj-guess repository, which you chose when you checked
out the source code via `svn checkout`.

```
cd /path/to/lpj-guess/build
cmake ..
cmake --build .
```

## Debugging

0. Open visual studio code
1. Open the repository if it is not already open ("File -> Open Folder")
2. Press F5 on the keyboard

This will prompt for an ozflux site to be simulated, and the name of an output
directory. All instruction files will be copied into this directory before the
simulation is run. All outputs will be written into this directory.

If using a keyboard without an F5 key, open the debug panel (view -> run, or
ctrl+shift+d, or cmd+shift+d on MacOS). Click on the green arrow near the top
of the debug panel (should be near the top-left of the screen, to the left of
some text which says "Run Site").

## Notes

### Linking against the NetCDF libs

If you haven't installed the netcdf libraries, you'll get a warning after running `cmake ..`:

```
find_package: not found: NETCDF_C_LIBRARY and NETCDF_INCLUDE_DIR. See LPJ-GUESS documentation.
NETCDF disabled.
```

This is a warning and not an error, so it's still possible to compile successfully despite this, it just means LPJ-Guess won't be able to read netcdf input files. This can be rectified by installing the netcdf package via homebrew.

If you encounter this warning on gadi, you'll need to load the netcdf module and tell cmake how to find it:

```bash
module load netcdf
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$NETCDF"
export ARCH=gadi # (recommended)
cmake ..
cmake --build .
```

## Troubleshooting

### Compiling/Debugging in vscode doesn't work

If attempting to run the model gives a message such as "select debugger", or
attempting to compile the model gives a message such as "select the build task
to run", then you have most likely opened the wrong directory (possibly a
subdirectory such as modules or framework?).

***When working with the model in visual studio code, always open the "main"
model directory via "File -> Open Folder".*** The "main" directory is the one
which contains other directories such as build, modules, framework, etc.

The name of the currently-opened directory is visible in the middle of the
window title bar at the top of the screen, next to the name of the currently
open file. When launched, vscode will usually open the previously open
workspace, and if that is the lpj-guess repository, there is no need to re-open
the same directory.
