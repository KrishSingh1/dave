Custom libraries used by LPJ-GUESS
==================================

LPJ-GUESS (versions from 2001-11-28) require the custom libraries gutil and
plib. These provide functions and classes of a purely technical nature which
serve to simplify the code of LPJ-GUESS itself.


Users building guess.exe or guess.dll under Windows
---------------------------------------------------

For Windows users, the libraries are provided in binary form. The binaries are
located under directory ./windows_bin/. Further information is given in
conjunction with instructions for building guess.exe or guess.dll, i.e. in the
files readme_windows.txt and readme_unix.txt in directory
../command_line_version and in the file ../windows_version/readme.txt.


Users building guess.exe under Unix
-----------------------------------

For Unix users, the libraries are provided as source code and MUST BE BUILT
BEFORE LPJ-GUESS IS BUILT. In general, however, users on the same system can
share a single copy of the library binaries.

If the libraries are not already available as binary archive (.a) files
somewhere on your system, do the following.

1. Create a directory where the libraries and include files will be stored.
   This directory may have any name, but the name /user/lib/ is assumed in
   these instructions.

2. Under /user/lib/ on your system, create a directory called include/

3. Copy the files under ./include/ (on this FTP area) to /user/lib/include/
   on your system.

4. Copy the following files from this FTP area to /user/lib/ on your system:
   - Makefile
   - the files under ./source/
 
5. You should now have the following or a similar hierarchy of files and
   directories on your system:

   /user/lib/
     gutil.cpp
     plib.cpp
     Makefile
     include/
       gutil.h
       plib.h
 
6. In a command shell, go to directory /user/lib/ on your system, then type
   "make"; this should build the binary archive files gutil.a and plib.a
   (two additional files, gutil.o and plib.o may also be created).

7. LPJ-GUESS may now be built as a Unix command-line executable by following
   the instructions in ../command_line_version/readme_unix.txt.


Ben Smith
2001-11-28
