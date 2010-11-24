===============================================================================
     Tools for building and running LPJ-GUESS for parallel runs on Simba
                    Portable Batch System (PBS) version
===============================================================================


                            PLEASE READ CAREFULLY


Adapting your job
-----------------

   These tools should work for most versions of LPJ-GUESS. However, you may
   have to make some changes to your current setup. The following requirements
   must be met EXACTLY. If not, the job will not work as expected or may not
   work at all.

   REQUIREMENT 1: GRIDLIST FILE
     Grid cells / site coordinates MUST BE read from one plain-text gridlist
     file (one row per grid cell). Other solutions, e.g. a lat-lon window,
     will not work with these tools. Your best bet is to modify your setup and
     implement a gridlist file.

   REQUIREMENT 2: GRIDLIST MUST END WITH NEWLINE
     Ensure that the last line of the gridlist file ends with a "newline"
     character, otherwise the last grid cell will not be simulated. If you are
     unsure whether this is the case, add a blank line at the end.

   REQUIREMENT 3: THERE MUST BE AN INS FILE

   REQUIREMENT 4: OUTPUT FILES MUST BE PLAIN TEXT FILES
     All output files must be plain text files with or without one leading row
     containing a header or column labels (or a blank first row). If your
     output files do not meet these conditions, they will still appear
     for subsets of gridcells in the runX subdirectory for each process
     (see Steps 8-9 below) but combined files containing output for all
     gridcells will be corrupted or may not be created.

   REQUIREMENT 5: FILE NAMES MUST BE SPECIFIED THE CORRECT WAY
     File names in the ins file and/or input/output module (guessio.cpp or
     equivalent) must be specified in the correct way: some files MUST be
     specified as a SIMPLE FILE NAME with NO DIRECTORY PART while others MUST
     be specified as a FULL PATH NAME INCLUDING A DIRECTORY PART.

     The rules are as follows. READ CAREFULLY AND FOLLOW EXACTLY.  

   - Full absolute pathnames must be given for all INPUT files the model
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
     specified, if required, in the submit script (see Step 8 below):

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
       param "file_cmass" (str "/scratch/fred/guessrun/cmass.out")
       file_flux="../guessrun/flux.out"



Building the model
------------------

1. You should be logged onto simba.nateko.lu.se.

2. Confirm that your code and setup meed the requirements specified above.
   PLEASE READ THESE CAREFULLY! If not, modify your setup so that they do.

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
   ftp> cd guess030124/parallel_version/simba/pbs
   ftp> get *
   ftp> exit

4. Copy all the LPJ-GUESS module and framework files into the build
   directory. However, do NOT copy over main.cpp (as it will write over the
   version you copied over in step 3) nor the libraries gutil.cpp or
   plib.cpp.

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

   cp /home/fred/guess_mpi/submit.sh .     [don't forget the '.' at the end!]

8. Open submit.sh in a text editor and set appropriate values for the
   following variables (NB: no space allowed after the "=" sign!):

                             -- READ CAREFULLY --
   
   NPROCESS = number of processes (nodes) you want for the parallel job.
              A typical number would be 4, the maximum that is sensible is
              8 (as there are 8 nodes on Simba). Nodes are reserved for
              exclusive use (i.e. all other jobs will be queued until the
              the current one is finished), so please be considerate to other
              users. Note that there is no point requesting more processes
              than there are grid cells in the gridlist file.

   WALLTIME = maximum wall (real) time for job in format hh:mm:ss
   
   INSFILE  = pathname (with or without directory part) to ins file FROM THE
              RUN DIRECTORY.

   GRIDLIST = pathname (with or without directory part) to gridlist file FROM
              THE RUN DIRECTORY.

   OUTFILES = list of LPJ-GUESS output files in single quotes, and separated
              by spaces (FILENAMES ONLY, INCLUDING EXTENSION, NO DIRECTORY).
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

10. The job number is displayed when you submit the job is Step 9. You can see
    all queued and running jobs using:

    qstat

    or (to see your own job(s) only)

    qstat -u <user-id>
    e.g. qstat -u fred

    Detailed status for a particular job is available using

    qstat -f <job-id>

    e.g. qstat -f 8410

11. A log file for each subprocess (guess.log) appears in its corresponding
    runX subdirectory. You can display the last few lines of each log file
    to check progress by entering the following command in the run directory:

    sh progress.sh

12. To cancel a running job, delete it using:

    qdel <job-id>
    e.g. qdel 8410

13. Final output appears in the run directory (see also Step 9 above).
    A log with some additional messages appears in stdout.log and in a file
    called something like guess.cmd.XXX in the run directory, where XXX
    is the job number from Step 9. In the event of errors these file may
    provide som clues.


Alternative options for advanced users
--------------------------------------

Janno has created a 'makescript' file that allows single- or dual-processor
nodes to be chosen in the submit.sh. This can be an advantage if you want to
ensure that all processes run equally fast (since the dual-processor nodes on
simba currently have faster clock speeds than the older single-processor nodes).
To use this feature, do the following in the BUILD directory:

     - delete the file makescript:

       rm makescript

     - rename makescript_flexproc to makescript:

       mv makescript_flexproc makescript

     - build a new submit script:

       make

     - repeat steps 6-9 above (see comments in submit.sh for new settings)


Troubleshooting
---------------

If it doesn't work:

  * Carefully reread "Adapting your job" above, especially REQUIREMENT 5.

  * Carefully reread step 8 under "Running the model" above.
  
  * Check for error messages in stdout.log and guess.cmd.eXXX (see Step
    13 above).

Specific problems and possible causes:

  PROBLEM: Error messages like "cannot remove XXX: No such file or directory"
    or "cannot create directory run1: File exists" when I run submit.sh
  POSSIBLE CAUSE: Probably there is no problem. If a job number followed by
    ".simba.nateko.lu.se" appears as the last line of output when you run
    submit.sh then the job has been successfully submitted to the batch queue.

  PROBLEM: Doesn't seem particularly fast for a supercomputer. I am trying
    this out for one grid cell.
  POSSIBLE CAUSE: There is no point running a parallel job for one grid cell,
    since the grid cells are divided equally between processes (nodes).
    See NPROCESS under Step 8 above.

  PROBLEM: Seems very slow, same grid cells appear several times in output.
  POSSIBLE CAUSE: Path name specified for gridlist file in ins file or i/o
    module. File name only (no directory part) must be specified.
    See REQUIREMENT 5 above.

  PROBLEM: Runs for several hours then stops before all grid cells have
    been processed.
  POSSIBLE CAUSE: Maximum time for job exceeded. You need to specify a
    longer wall clock time in the submit script. See Step 8 above.

  PROBLEM: Job seems to run and finish normally, but no output files appear
    in the run directory
  POSSIBLE CAUSE: Incorrect or missing output file list in submit script.
    Should be a space-delimited list of file names including extensions
    and with no directory part. See Step 8 above.

  PROBLEM: Output files are unreadable or corrupted
  POSSIBLE CAUSE: The model (input/output module) is producing output files
    in the wrong format. Output files must be plain text files with or without
    ONE initial row containing header information, column labels or a blank
    line. See REQUIREMENT 4 above.

  PROBLEM: Output is missing for the last grid cell in the gridlist file
  POSSIBLE CAUSE: The last line of the gridlist file does not end with a
    newline character (see REQUIREMENT 2 above). Add a blank line at the end
    of the gridlist file.

  PROBLEM: Job seems to enter the batch queue but crashes immediately after
    start.
  POSSIBLE CAUSE: Format errors or incorrect file names specified in
    submit.sh. (NB: no space allowed after "=" in variable assignments!)
    See Step 8 above.

  PROBLEM: Job appears in batch queue but never seems to start
  POSSIBLE CAUSES: You are requesting more processes than there are nodes
    available. Your job will start when there is one free node for each
    process. Simba has (at the time of writing) 8 nodes, but some may be
    reserved/in use by other users. Try requesting fewer processes for your
    job. See Step 8 above.


Ben Smith
2005-03-24
