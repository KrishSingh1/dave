//////////////////////////////////////////////////////////////////////////////
// Building and setting up LPJ-GUESS for parallel runs on Lunarc/Milleotto

Building an executable:

In addition to framework and module source code (.cpp) and header (.h) file,
the following files in this directory are required:

  - main.cpp
  - append.cpp
  - splitgrid.cpp
  - makescript
  - Makefile

(all except append.cpp and splitgrid.cpp are specific to Milleotto and will
not work on other Lunarc systems such as docenten).

* Copy all the above files, including framework and module source and header
  files to a directory under /disk/global/(your user id). This is called the
  BUILD DIRECTORY

* Add the following declaration near the top of guess.h (after the
  #include's):

  void dprintf_mon(xtring format,...);

* Add the following statement in function getstand() in guessio.cpp:
  (required to provide limited 'log' output while GUESS is running - 
  output from ordinary dprintf statements will not be available until the
  end of the run): 

  dprintf_mon("\nCommencing simulation for stand at (%g,%g)\n",
     gridlist.getobj().lon,gridlist.getobj().lat);

* Make any necessary amendments to the Makefile (e.g. non-standard filenames)

* Load some modules required for the build:

  module load pgi
  module load mpich-pgi6

* Build the model and supplementary programs and scripts (append, splitgrid,
  submit.sh):

  make


Running the model in parallel mode:

* Create a directory under /disk/global/(your user id) to which output from
  the model run should be written. This is called the RUN DIRECTORY.

* Go to the run directory, e.g.

  cd run

* Copy the file 'submit.sh' (created during the build above) from the build
  directory to the run directory.

* Make necessary amendments to settings listed and explained near the top of
  submit.sh

* Submit the parallel job:

  sh submit.sh

* Job status may be checked by:

  qstat -u (your user id)

* Progress (grid cells started) may be checked while the job is running by
  going to the run directory and typing:

  sh progress.sh

* If the run is successful, output is written to the run directory


Further information:

More general information is available in ../docenten/readme.txt


Ben Smith
2007-08-16
