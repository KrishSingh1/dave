#!/bin/bash

# Delta (deltareports)
# A script generating delta reports: absolute and relative differences between the LPJ-GUESS traditional benchmarks
# and a the correspoding ones of a refernce model run.
# The script takes 2 arguments: a label, and a path, for the folder containing the benchmarks to compare with.


resultfile="deltareports.vs-${1}.result"

echo "Summary-processing: Deltareports... Arguments = $@"
echo "$0"
echo "The summarize script's pwd -P = $(pwd -P)"
cd ..
#which deltareports | tee -a delta/$resultfile
deltareports -c $@ | tee -a delta/$resultfile


