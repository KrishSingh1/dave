LPJ-GUESS Version 2001-11-28 : Setup as a Windows command line executable
=========================================================================

This file explains how to set up LPJ-GUESS as a command line executable
(console application) to run in a command shell (MS-DOS window) under Windows.
These instructions assume you have Microsoft Visual C++ / Developer Studio
Version 6.

1. In Visual C++/Developer Studio, create a new workspace called Guess. This
   will build LPJ-GUESS as a command line executable.

    --> Choose File|New. On the "Project" tab

 	--> choose "Win32 console application"
        --> set "Project name" and "Location". You may use any project name
            you wish. In these instructions, the project name Guess, and
            workspace directory C:\Guess\, are assumed
        --> press "OK"
        --> press "Finish"

2. Transfer the following files by FTP to the Guess workspace directory:
   (NB: ASCII mode):

   - all the .cpp and .h files under ../modules/
   - the .cpp and .h file under ../framework/
   - the file main.cpp in the current directory (../command_line_version/)

3. Create a directory on your PC that will store the binary libraries and
   library include files required by LPJ-GUESS. In these instructions, this
   directory is assumed to be C:\lib\, but you may choose any location you
   wish. Under this directory create an additional directory called "include".

   -->  Transfer (in ASCII mode) the files under ../libraries/include/
        to the directory C:\lib\include\
   -->  Transfer (in binary mode) the .lib files under
        ../libraries/windows_bin/ to the directory C:\lib\

4. Return to Visual C++, reopen workspace "Guess" if necessary

   -->  Using the menu option Project|Add to project...|Files, add the
        following files to the Guess project:

        - all .cpp and .h files in the workspace directory (i.e. the module
          and framework files, plus main.cpp)
        - the library files gutil.lib and plib.lib, located in C:\lib\

   -->  Choose Tools|Options, go to the "Directories" tab. In the window
        labelled "Show directories for" choose "Include files", then add the
        pathname of the directory containing the library include files, i.e.
        C:\lib\include

   -->  Choose Build|Build guess.exe. This should build the binary executable
        file in the Debug subdirectory. The compiler and/or linker may report
        warnings, but there should be no errors

   -->  If GUESS was successfully built, you should now be able to run it at a
        command prompt. You can also run GUESS from within Visual C++/
        Developer Studio, by choosing the menu option Build|Execute. In both
        cases you need to specify the pathname of an instruction script (ins)
        file as a command line argument. From Developer Studio, you specify
        the argument as part of the project settings:
        
        --> choose Project|Settings
        --> select the "Debug" tab
        --> under "Program arguments" enter the name or pathname of an ins
            file

        NOTE: the project settings have no effect when guess.exe is run
              directly at a command prompt. In this case, the ins file should
              be specified as the (only) command line argument; e.g.

              guess demo.ins

        You may also run guess with "-help" as the command line argument,
        instead of the name of an ins file. This results in output of a list
        of keywords recognised in the ins file.

The default input/output module (guessio.cpp) distributed with LPJ-GUESS
requires Cramer & Leemans (unpublished) temperature, precipitation and
sunshine files, the LPJ soil code file, a file containing a list of grid cells
to simulate, and an appropriate ins file. These files can be found under
../data.

Ben Smith
2001-11-28
