===============================================================================
     Tools for building and running LPJ-GUESS for parallel runs on LUNARC
===============================================================================


                           PLEASE READ CAREFULLY

 
IMPORTANT: The tools provided in this directory and the method described in
  this file will work only on the 'docenten' system and not on the system
  'toto7/whenim64'. This applies both to building and running the model.


What sort of jobs will this work for?
-------------------------------------

   It should work for most versions of LPJ-GUESS, provided that:

   - Grid cells / site coordinates are read from one plain-text gridlist file
     (one row per grid cell). Other solutions, e.g. a lat-lon window, will
     not work with these tools. Your best bet is to modify your setup and
     implement a gridlist file. IMPORTANT: ensure that the last line of the
     gridlist file ends with a "newline" character, otherwise the last grid
     cell will not be simulated. If you are unsure whether this is the case,
     add a blank line at the end.

   - There is an ins file

   - Full absolute pathnames are given for all INPUT files the model
     requires (NB: except the gridlist file, see next point), whether these are
     specified in the ins file or in the code:

     OK:
       param "file_cru" (str "/home/ben/archive/env/cru/cru_1901-1998.bin")
       file_soil="/home/ben/archive/env/cralee/soils_lpj.dat";

     NOT OK:
       param "file_cru" (str "cru_1901-1998.bin")
       file_soil="../guess/soils_lpj.dat";

   - The gridlist file name must be given as a SIMPLE FILE NAME with NO 
     DIRECTORY PART, whether it is specified in the ins file or in the code.
     However, a full or relative pathname for the gridlist file may be
     specified, if required, in the submit script (see Step 9 below):

     OK:
       param "file_gridlist" (str "gridlist.txt")
       file_gridlist="gridlist.txt";

     NOT OK:
       param "file_gridlist" (str "/home/ben/guess/gridlist.txt")
       file_gridlist="../gridlist.txt";

   - OUTPUT file names are specified as FILE NAMES ONLY with NO DIRECTORY
     PART, whether these are specified in the ins file or in the code:

     OK:
       param "file_cmass" (str "cmass.out")
       file_flux="flux.out";

     NOT OK:
       param "file_cmass" (str "/data/global1/fred/guessrun/cmass.out")
       file_flux="../guessrun/flux.out"

   - All output files are plain text files with or without one leading row
     containing a header or column labels (or a blank first row). If your
     output files do not meet these conditions, they will still appear
     for subsets of gridcells in the runX subdirectory for each process
     (see Step 10 below) but combined files containing output for all
     gridcells will be corrupted or may not be created.


Building the model
------------------

1. Confirm that you are logged onto docenten.lunarc.lu.se. If not, log out
   and start again.

2. Confirm that your code and setup meed the requirements specified above.
   READ THESE CAREFULLY! If not, modify your setup so that they do.

3. Load the pgi (compiler) and mpi (message pasting interface) modules:

   module add pgi
   module add mpich-pgi5

4. Create a directory in which to build the model. Copy all the files in
   the present directory (the one containing this readme.txt file) to the
   new directory you created. This will be called the BUILD DIRECTORY. In
   these instructions the build directory will be assumed to be called
   /home/fred/guess_mpi:

   cd /home/fred/
   mkdir guess_mpi
   cd guess_mpi
   sftp fred@stormbringerii.nateko.lu.se
   ftp> cd /user/home/ben/public_html
   ftp> cd guess030124/parallel_version/lunarc/docenten
   ftp> get *
   ftp> exit

5. Copy all the LPJ-GUESS module and framework files into the build
   directory. However, do NOT copy over main.cpp (as it will write over the
   version you copied over in step 4) nor the libraries gutil.cpp or
   plib.cpp.

   cp /someplace/guess.cpp . 
   cp /someplace/canexch.cpp .
   (etc)
   
   cp /someplace/guess.h .
   cp /someplace/canexch.h .
   (etc)

6. Build the model and portable submit script:

   make


Running the model
-----------------

7. Create or go to a directory in which to run the model. This will be
   called the RUN DIRECTORY. It should preferably be somewhere in the
   scratch area under /disk/global1/, NOT under /home. There is a 
   workspace with your user-id under /disk/global1:

   cd /disk/global1/fred
   mkdir guessrun
   cd guessrun

8. In step 6 you created submit.sh, a portable bash shell for running the
   model as a parallel batch job. Copy submit.sh to the above directory:

   cp /home/fred/guess_mpi/submit.sh .           [don't forget the '.'!!]

9. Open submit.sh in a text editor and set appropriate values for the
   following variables:

   NPROCESS = number of processes (nodes) you want for the parallel job.
              A typical number would be 8, the maximum permitted is 32.
              The more nodes, the faster the job will run once it starts.
              However if you request many nodes your job may wait a long time
              for these to become available. Note that there is no point
              requesting more processes than the number of grid cells in the
              grid list file.

   WALLTIME = maximum wall (real) time for job hh:mm:ss. This must be long
              enough for your job to run to completion, but the shorter it is,
              the sooner your job will start.

   EMAIL    = to which notifications of job completion/failure should be sent.
              Please be sure to enter YOUR e-mail address here.

   INSFILE  = pathname (with or without directory part) to ins file from the
              run directory.

   GRIDLIST = pathname (with or without directory part) to gridlist file from
              the run directory.

   OUTFILES = list of LPJ-GUESS output files in single quotes, and separated
              by spaces (filenames only, including extension, no directory).
              Only plain text files should be included in this list. If you
              insist on producing other kinds of output (e.g. binary files),
              they will appear (for a subset of grid cells) in the runX
              subdirectories, but should not be included in this list, and
              will not be appended together at the end of the run.

10. Submit the job to the batch queue:

   sh submit.sh

   Once the job starts it will create a number of subdirectories called runX
   (X=process number), one for each process/node in the parallel job. A
   roughly equal chunk of the gridlist file and a copy of the ins file appears
   in each subdirectory. Output appears at the end of the run only in each
   subdirectory (for the gridcells listed in the gridlist file there) and
   in the run directory (all gridcells).


Checking progress etc.
----------------------

11. You can see your job(s) in the batch queue using:

    qstat -u <your-user-id>
    e.g. qstat -u fred

12. Once the job is running you can upload the current guess.log for each
    process to its corresponding runX directory by running the progress.sh
    shell in the directory from which you submitted the job. The last few
    lines of each log file are also displayed on the screen:

    sh progress.sh

13. To cancel a waiting or running batch job, find out its id number using
    qstat (Step 11 above), then delete it using:

    qdel <job-id>
    e.g. qdel 8410

14. Final output appears in the run directory (see also Step 10 above). 


Ben Smith
2005-03-24

