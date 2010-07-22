LPJ-GUESS setup as a command line executable on Simba
=====================================================

This version will produce graphical output from “plot()” function calls in the
model code if an X11 graphics server (X-Server) is available and if graphical
output is not disabled by the “-nogra” option (see below).

1. Create a directory to store the source code and header files, make file and
   (once built) the binary executable. This is referred to in these
   instructions as the build directory.

2. Transfer the following files to the build directory

   - all the .cpp and .h files under ../modules/
   - the .cpp and .h file under ../framework/
   - the files main.cpp and Makefile in the current directory
     (../command_line_version/simba)

3. The make file (Makefile) contains information about pathnames to libraries
   and include files which are required to compile LPJ-GUESS. Ensure that the
   following files are available in the specified directories; change the
   relevant path names in the make file if necessary:

   - the library (binary object) files gutil.o and plib.o
   - the include files gutil.h and plib.h

4. Build the executable file (guess):

   - go to the build directory
   - at the command prompt, type: make

5. If guess was successfully built, you should be able to run it at a command
   prompt. The syntax for this is:

   guess [options] <instruction-file>

   Options: -log <file>       specifies log file name
            -nolog            no log file
            -display <dest>   address or ip-number to X windows server
            -mono             monochrome graphical output
            -width <pixels>   graph window width
            -height <pixels>  graph window height
            -nogra            suppress graphical output
            -wait             wait for keyboard at end of run (to preserve graphics)
            -help             list ins file parameters  guess lpj.ins

   For example:
   
   guess -nogra lpj_030124.ins

   Demonstation instruction files and the input data they require may be found
   under ../data.

   
Ben Smith
2005-01-26
