#!/bin/bash

# Delta (deltareports)
# A script generating delta reports: absolute and relative differences between the LPJ-GUESS traditional benchmarks
# and a the correspoding ones of a refernce model run.
# The script takes 2 arguments: a label, and a path, for the folder containing the benchmarks to compare with.


resultfile="deltareports.vs-${1}.result"

echo "Exo-processing: Deltareports... Arguments = $@"
deltareports $@ | tee -a $resultfile


############## Old

# Decision by Johan and Stefan 2022-02-11
# to rename postpostprocess-files and folders to exoprocess.
# because postpostprocess is impossible to pronounce.
# exo meaning "on the outside" etc.
# e.g.:
#      benchmarks\postpost\delta\postpostprocess.sh"
#  --> benchmarks\exoprocesses\delta\exoprocess.sh"
#                  ^                  ^

#echo "Delta tool" | tee postpostprocess.dummy-output
#echo "Arguments = $@" | tee -a postpostprocess.dummy-output
#echo "Argument 1 = $1" | tee -a postpostprocess.dummy-output
#echo "Argument 2 = $2" | tee -a postpostprocess.dummy-output
#echo "pwd -P = $(pwd -P)" | tee -a postpostprocess.dummy-output
