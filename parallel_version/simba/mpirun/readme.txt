===============================================================================
     Tools for building and running LPJ-GUESS for parallel runs on Simba
                         Non-queue (mpirun) version
===============================================================================


*******************************************************************************
*                                                                             *
*             STOP GO BACK! THIS VERSION NOT TO BE USED ON SIMBA              *
*                                                                             *
*    As of 2005-03-16 the Portable Batch System (PBS) has been                *
*    installed on Simba. All parallel runs should now be performed using PBS  *
*    rather than directly via mpirun (i.e. the tools in this directory).      *
*                                                                             * 
*    Instructions for building and running LPJ-GUESS on Simba using PBS are   *
*    available in directory ../pbs                                            *
*                                                                             *
*    Ben Smith, 2003-03-16                                                    *
*                                                                             *
*******************************************************************************



                           PLEASE READ CAREFULLY


What sort of jobs will this work for?
-------------------------------------

   It should work for most versions of LPJ-GUESS, provided that:

   - Grid cells / site coordinates are read from one plain-text gridlist file
     (one row per grid cell). Other solutions, e.g. a lat-lon window, will
     not work with these tools. Your best bet is to modify your setup and
     implement a gridlist file.

   - There is an ins file

   - Full absolute pathnames are given for all INPUT files the model
     requires (NB: except the grid list file), whether these are specified in
     the ins file or in the code:

     OK:
       param "file_cru" (str "/home/ben/archive/env/cru/cru_1901-1998.bin")
       file_soil="/home/ben/archive/env/cralee/soils_lpj.dat";

     NOT OK:
       param "file_cru" (str "cru_1901-1998.bin")
       file_soil="../../archive/env/cralee/soils_lpj.dat";

   - OUTPUT file names are specified as file names only with no directory
     part, whether these are specified in the ins file or in the code:

     OK:
       param "file_cmass" (str "cmass.out")
       file_flux="flux.out";

     NOT OK:
       param "file_cmass" (str "/scratch/fred/guessrun/cmass.out")
       file_flux="../guessrun/flux.out"

   - All output files are plain text files with or without one leading row
     containing a header or column labels (or a blank first row). If your
     output files do not meet these conditions, they will still appear
     for subsets of gridcells in the runX subdirectory for each process
     (see Steps 8-9 below) but combined files containing output for all
     gridcells will be corrupted or may not be created.


Building the model
------------------

1. You should be logged onto simba.nateko.lu.se.

2. Confirm that your code and setup meed the requirements specified above.
   If not, modify your setup so that they do.

3. Create a directory in which to build the model. Copy all the files in
   the present directory (the one containing this readme.txt file) to the
   new directory you created. This will be called the BUILD DIRECTORY. In
   these instructions the build directory will be assumed to be called
   /home/fred/guess_mpi:

   cd /home/fred/
   mkdir guess_mpi
   cd guess_mpi
   sftp fred@stormbringerii.nateko.lu.se
   ftp> cd /user/home/ben/public_html
   ftp> cd guess030124/parallel_version/simba
   ftp> get *
   ftp> exit

4. Copy all the LPJ-GUESS module and framework files into the build
   directory. However, do NOT copy over main.cpp nor the libraries
   gutil.cpp or plib.cpp.

   cp /someplace/guess.cpp . 
   cp /someplace/canexch.cpp .
   (etc)
   
   cp /someplace/guess.h .
   cp /someplace/canexch.h .
   (etc)

5. Build the model and portable submit script:

   make


Running the model
-----------------

6. Create or go to a directory in which to run the model. This will be
   called the RUN DIRECTORY. It should preferably be somewhere in the
   scratch area under /scratch, NOT under /home. There is a workspace with
   your user-id under /scratch:

   cd /scratch/fred
   mkdir guessrun
   cd guessrun

7. In step 5 you created submit.sh, a portable bash shell for running the
   model as a parallel batch job. Copy submit.sh to the run directory:

   cp /home/fred/guess_mpi/submit.sh .          [don't forget the '.' !!]

8. Open submit.sh in a text editor and set appropriate values for the
   following variables:

   NPROCESS = number of processes (nodes) you want for the parallel job.
              A typical number would be 4, the maximum that is sensible is
              8 (as there are 8 nodes on Simba). More processes will improve
              run time if there are resources available, but if demand from
              other users and jobs is high, it will just slow everybody down.

   INSFILE  = pathname to ins file from the run directory.

   GRIDLIST = pathname to gridlist file from output directory.

   OUTFILES = list of LPJ-GUESS output files in single quotes, and separated
              by spaces (filenames only, including extension, no directory).
              Only plain text files should be included in this list. If you
              insist on producing other kinds of output (e.g. binary files),
              they will appear (for a subset of grid cells) in the runX
              subdirectories, but should not be included in this list, and
              will not be appended together at the end of the run.

9. Start the job:

   sh submit.sh

   Once the job starts it will create a number of subdirectories called runX
   (X=process number), one for each process/node in the parallel job. A
   roughly equal chunk of the gridlist file and a copy of the ins file appears
   in each subdirectory. Output appears in each subdirectory (for the
   gridcells listed in the gridlist file there) and, at the end of the run,
   in the run directory (all gridcells).


Checking progress etc.
----------------------

10. The process number is displayed when you submit the job is Step 9. The
    process name is "mpirun". You can see all processes under your user-id
    using:

    ps -u <your-user-id>
    e.g. ps -u fred

11. A log file for each subprocess (guess.log) appears in its corresponding
    runX subdirectory. You can display the last few lines of each log file
    to check progress by entering the following command in the run directory:

    sh progress.sh

12. To cancel a running job, delete it using:

    kill <process-id>
    e.g. qdel 8410

    You can kill all processes under your name on the system using the
    following, which will also log you out:

    kill -9 -1

13. Final output appears in the run directory (see also Step 9 above).
    A log with some additional messages appears in stdout.log in the
    run directory. In the event of errors this file may provide som clues.


Ben Smith
2005-01-26
