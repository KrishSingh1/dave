The tool scripts here can be submitted to the HPC cluster nodes by the script benchmarks.
Run ./benchmarks for help.
But for help on the individual scripts here, please refer to the commenting at the 
top of each script: LPJ-GUESS code directory/benchmarks/summaries/<tool>/summarize.sh

Note that the benchmarks functionality of LPJ-GUESS only works on linux.

Summary jobs is a mechanism within the benchmark command to let specialized jobs 
start *after* all postprocessing jobs are finished.
Therefore the summary jobs can be used to summarize the result of several 
or all of the postprocess results of the benchmark run.
Several summary jobs can be started within one overall benchmark command run.
The summary jobs run in parallel independent of each other.

So for example:
The benchmark run can be run with all benchmark tests...
 ./benchmarks <output dir>
... or run with only a few benchmark tests:
 ./benchmarks -i "global crop_global" <output dir>
But since the summary jobs run in parallel independently of which  
benchmark model runs that has been run (and their associated postprocessings),
they need to be individually specified in the benchmarks command call:
 ./benchmarks -s "diff <label> <comparison output>" <output dir>

In the example above -s is the command option that tells the benchmark command
to include a summary job in the run (i.e. in is submitted jobs to the cluster),
and the summary job to include in this example is the "diff" summary command.
The "diff" summary command, takes 2 arguments: a label and a path to the reference 
output to diff the current run's output with (e.g. a previous output from trunk).
The purpose of the label is to distinguish different diff summary jobs' results
from each other (e.g. comparison with trunk_r13354 in one summary job and 
european_applications_rxxxxx in another summary job).
Note that the summary command and its arguments are enclosed within quotes after 
the -s option switch.

Several summary jobs can be run in the same benchmarks run:
 ./benchmarks -s "diff <arguments...>" -s "tellme <arguments...>" <output dir>

The results of the summary jobs are written to the <output dir> of the benchmarks 
command, that is, in the same folder as the crop_global and global (etc) benchmarks'
results are written to.
E.g.:
ls <output dir>:
crop_global diff global tellme
Thus the results of a summary job will be written to a subfolder named as the summary tool.
The summary process will run inside that folder. I.e. slurm logs and error logs for the 
summary job will also be found withing that folder.
 
The currently available summary jobs.
All summary job scripts are located in the LPJ-GUESS code directory
benchmarks/summaries/
Specific info for each summary command can be read in a comment in the top of 
each script: LPJ-GUESS code directory/benchmarks/summaries/<tool>/summarize.sh
The currently available summary tools are:
* logreport: summarizes the guess.log log files of all the model runs of
  the benchmark run. Standard output is removed so that only unexpected output 
  such as warnings stand out. Expected files that are missing are also higlighted.
  -s "logreport <reference output>"
* diff: diffs the .out files and tslices of all the model runs of
  the benchmark run with some reference output (e.g. a trunk model run).
  Usefull if you expect that your changes will not affect model output and wnat to
  confirm that quickly.
  -s "diff <reference label> <reference output dir>"
* delta: like diff above, but used instead when you *do* expect changes in model output
  and want to see how large a difference (delta) your model changes has made to the output.
  This summary tool generates one delta-report inside each benchmark's outputfolder.  
  -s "delta <reference label> <reference output dir>"
* tellme: The tellme summary tool replaces the traditional benchmarks reports.
  It is a summary job because it reads output from several benchmark tests to 
  make one report of them all. The resulting report is written to the <output dir>/tellme
  folder.
  -s "tellme <this-job-label> <reference label> <reference output dir>"
Planned, not implemented yet:
* std-summaries: a summary meta-job that runs all 4 summary jobs above, so that you 
  don't have to specify <ref-label> <reference output dir> many times in your
  benchmarks command.
  -s "std-summaries <thisjob-label> <reference-label> <reference output dir>"
  or, if you'd like to exclude the heavy tellme summary job (i.e. only the 3 first jobs):
  -s "std-summaries <ref-label> <reference output dir>"

How to write a new summary tool:
Easy! 
* Make a new folder in LPJ-GUESS code directory/benchmarks/summaries/<tool>/
  in your branch. Name the new folder to the name of your new sumary tool. let's say "NEW".
* In the new folder NEW make file summarize.sh
  Easiest that you copy the summarize.sh from one of the other tools, e.g. the logreport one.
* Edit the comment section (this is the help text)
  and update the last part that calls the actual work script.
  The currently existing summary tools consists of 
  - a summarize.sh script that is called by the benchmarks command
  - and a work script that does the actual work. It is called by summarize.sh.
  You can do like that or have all you code in the summarize.sh script.
* Important: the pwd (current dir) of execution is down in the NEW folder.
  So your new script will need to look for the benchmark runs output one level up, in NEW/../
  But you will will want your summary tool's output to be written in NEW/
  Look in e.g. tool "diff" for how that was solved there.

Planned additions:
* Mentioned above, the std-summaries summary tool.
* A help command option in the benchmarks command that prints the top comment of each 
  summary tool.
* End jobs.  Can be used e.g. to zip all outfiles at the end of a benchmark run when all
  postprocessing and all summary jobs have finished.
  Summary jobs is a mechanism within the benchmark command to let specialized jobs 
  start after all *postprocessing jobs* are finished.
  In contrast, end jobs is a mechanism to let specialized jobs start after 
  *all* jobs are finished, not only postprocessing jobs.
  I.e. each end jobs is run one at a time, not starting until all other jobs,
  including a previous end job has finished.
