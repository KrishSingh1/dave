The tool scripts here can be submitted to the HPC cluster nodes by the script benchmarks.
run ./benchmarks for help.
But for help on the individual scripts here, please refer to the commenting at the top of each script.

For each script specified by the user, benchmarks will create a subfolder named as the tool script.
And the script process will run inside that folder.
So to access the model benchmarks run's output, the tool script needs to use the relative path ../
in order to access to top level of the benchmarks output folder.

Note that if you write a new script, say NEW/:
* if you do cd .. in your script in NEW/,then, if you redirect output to a file, 
  you will need to redirect it to >NEW/result.file
* Write a comment beginning on line 3 that describes your script.
