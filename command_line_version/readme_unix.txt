LPJ-GUESS Version 2001-11-28 : Setup as a Unix command line executable
=========================================================================

This file explains how to set up LPJ-GUESS as a command line executable to run
under Unix. These instructions assume you have a C++ compiler called CC with
standard compiler options.

1. Make sure the custom libraries gutil and plib are available on your system
   as binary archive files (gutil.a and plib.a). If not, you will have to
   build them first (see ../libraries/readme.txt)

2. Create a directory to store the source code and header files, make file and
   (once built) the binary executable. This is referred to in these
   instructions as the workspace directory.

3. Transfer the following files to the workspace directory

   - all the .cpp and .h files under ../modules/
   - the .cpp and .h file under ../framework/
   - the files main.cpp and Makefile in the current directory
     (../command_line_version/)

4. The make file (Makefile) contains information about pathnames to libraries
   and include files which are required to compile LPJ-GUESS. Ensure that the
   following files are available in the specified directories; change the
   relevant path names in the make file if necessary:

   - the library (binary archive) files gutil.a and plib.a
   - the include files gutil.h and plib.h

5. Build the executable file (guess):

   - In a unix command shell, go to the workspace directory
   - At the command prompt, type: make

6. If guess was successfully built, you should be able to run it at a command
   prompt. An instruction script (ins) file name must be specified as a
   command-line argument, for example:

   guess lpj.ins

   You may also run guess with the -help option to obtain a list of keywords
   recognised in the ins file:

   guess -help

The default input/output module (guessio.cpp) distributed with LPJ-GUESS
requires Cramer & Leemans (unpublished) temperature, precipitation and
sunshine files, the LPJ soil code file, a file containing a list of grid cells
to simulate, and an appropriate ins file. These files can be found under
../data.

Ben Smith
2001-11-28
