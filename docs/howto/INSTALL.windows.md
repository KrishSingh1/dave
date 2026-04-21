# LPJ-GUESS DAVE Windows Setup Guide

## Preamble

This guide explains how to compile, run, and debug LPJ-GUESS DAVE on windows,
using the msys2 toolchain + vscode. It's possible to instead use Microsoft
compilers (MSVC) and visual studio, and while that will work, it's beyond the
scope of this guide (ie: you're on your own).

## Software Installation

- Download and install filezilla
- Download and install sourcetree (skip bitbucket config, disable mercurial)
- Download and install visual studio code (not visual studio)
- Download and install msys2 (I recommend using the default installation path)

## Software Configuration/Setup

### MSYS2

Run the following command from inside msys2:

```bash
pacman -S mingw-w64-x86_64-{gcc,cmake,netcdf,gdb,make,ninja}
```

### Update PATH to include MinGW64

Go to the start menu and search for "environment", and click on "Edit
environment variables for the current user". Double click on "Path", then click
"New". Here you will need to use the path into which you installed msys2. By
default this was `C:\msys64`, but if you installed it in a different location,
use that location instead.

Add the following path: `C:\path\to\msys64\mingw64\bin`.

E.g. for the default installation options: `C:\msys64\mingw64\bin`.

Click OK to save the changes, and click OK on the environment variables window
as well.

To make sure this worked, open a **new** command prompt and run `gcc`. If you
get a "command not found" error, then something has gone wrong - talk to Drew.
If you get something like "fatal error: no input files", then you're all set.

If Visual Studio Code is installed and running, close it and reopen it.

### SSH Key Generation

**DO NOT do this if you already have a key.**

To check if you have a key, navigate to `C:\Users\username\.ssh`. If this
directory doesn't exist, you don't have a key. If you have a file in here called
something like `id_ed25519`, `id_rsa`, etc, then you already have an SSH key.

**Never** overwrite this key unless you know what you're doing.

If you don't have an SSH key, open a command prompt and run the following
command:

```bash
ssh-keygen -t ed25519
```

Just click "enter" at each prompt to generate a key with the default name and
without a passphrase.

Your key file comes in two parts, stored as two separate files in the
`C:\Users\username\.ssh` directory:

- `id_ed25519`                       (*Private key*)
- `id_ed25519.pub`                   (*Public key*)

In windows, it can be difficult to distinguish between the two files, so it can
be useful to enable the display of file path extensions in windows explorer.

### GitLab Account

If you don't already have a gitlab account, navigate to the
[Lund Gitlab](https://stormbringer4.nateko.lu.se/) in a browser and create one.

### GitLab SSH Key

Go to your
[preferences page](https://stormbringer4.nateko.lu.se/-/profile/preferences) in
gitlab, and click "SSH Keys" on the left, click "Add new key".

Copy the contents of your public key (`C:\Users\username\.ssh\id_ed25519.pub`)
and add it to your GitLab account. Make sure you copy the file ending with .pub,
which is your public key, as opposed to the file with the same name but no file
extension, which is your private key and should be treated like a password (ie
not shared with anyone).

To get the contents of your public key, you can open it in a text editor such as
notepad, or print it from a command prompt (e.g.
`type C:\Users\username\.ssh\id_ed25519.pub`)

### FileZilla

Go to the site manager (top-left, beneath the "file" button). Click "New site"
to register a new server with the following details:

- Protocol: SFTP - SSH File Transfer Protocol
- Host: hie-storage.sstars.ws
- Port: leave empty
- Logon Type: Key file
- User: Your hie-storage username (e.g. u90123456)
- Key file: Your **private** key (e.g. `C:\Users\username\.ssh\id_ed25519`)

Click OK to save this server configuration. To later connect to hie-storage,
select this server from the drop-down list to the right of the "site manager"
button.

### Visual Studio Code

Open Visual Studio Code, hit ctrl+shift+x to open the extension manager, and
install the following extensions:

- C++
- Native Debug

## Obtaining Source Code

Open SourceTree and click "clone". Use the following details:

Remote URL: `git@stormbringer4.nateko.lu.se:lpj-guess-developers/LPJ-GUESS.git`

Destination Path: `C:\some\path\on\your\computer`

The code will be downloaded into the destination path, so choose wherever works
best for you.

## Obtaining Climate Data

Climate data must be downloaded from hie-storage. If you don't have an account,
talk to Drew.

On hie-storage, the following directory must be downloaded:

```
/data/hiestorage/WorkingData/MEDLYN_GROUP/PROJECTS/dynamics_simulations/lpj-guess/forcing/ozflux/processed
```

This directory must be downloaded into your source code directory, inside
`<source code directory>\benchmarks\ozflux`, and it must be renamed to
forcing. So you should have a directory named
`<source code directory>\benchmarks\ozflux\forcing\met`.

## Running the Model

Open Visual Studio Code. If it was open before you edited the environment
variables above, close it and open it again.

Open the "Run and Debug" tab in the sidebar on the left. The icon is an
equilateral triangle with a bug on it. This can also be found by using the
keyboard shortcut "ctrl+shift+D".

There should be a dropdown which says "run site". This drop-down list contains
all of the ozflux sites at which the model is configured to run by default. To
run a particular site, select it from the dropdown list and hit the green play
button to the left of the dropdown list.

The model will now begin to run. You can watch the progress in the terminal
window at the bottom of the screen.

After the model has run, the outputs will be stored in:

`<source code directory>\benchmarks\ozflux\<site name>\out`

Running the model like this will automatically compile the code if necessary.

## Configuring the Model

This section could be arbitrarily long, so here I will briefly describe the
input file tree and list the main areas of interest.

For navigating these files, I recommend using the Visual Studio Code builtin
file explorer (top button in the left-hand sidebar, above run and debug, or use
the keyboard shortcut "ctrl+shift+e"). Alternatively, if you want to go to a
specific file, use the keyboard shortcut "ctrl+p", type the file name, and hit
enter.

All of the input files related to the ozflux sites are in the
`benchmarks\ozflux` directory. In here is one directory for each site. Each
site-level directory contains its own instruction file, named `site.ins`. These
site-level .ins files are very short but import other instruction files from
the ozflux directory, which contain the parameters common to all sites:

`settings.ins`: Contains most top-level switches, and, near the bottom of the
file, controls at a high level which kinds of PFTs are included in the
simulation.

`outputs.ins`: Contains all available output files which can be generated by the
model, most of which are disabled by default. To enable or disable an output,
simply comment/uncomment it. The comment character is an exclamation mark (!).

`trees.ins`: Contains parameters for all tree PFTs, as well as fine-grained
control over which specific trees are included in the simulation.

`c3g.ins`: Contains parameters for all C3 grass PFTs.

`c4g.ins`: Contains parameters for all C4 grass PFTs.

`pft_common.ins`: Contains parameters common to all PFTs.
