#!/bin/bash

# Delta (deltareports.sh -> deltareport -> deltareport_sub_imgs.sh)
# A script generating delta reports: absolute and relative differences between the LPJ-GUESS traditional benchmarks
# and a the correspoding ones of a refernce model run.
# The script takes 2 arguments: a label, and a path, for the folder containing the benchmarks to compare with.

toolscript="deltareports.sh"
resultfile="deltareports-vs-${1}.result"

if [ $# -ne 2 ]; then
  echo "Error: the summaries/delta/summarize.sh script must be called with exctly 2 arguments."
  echo "The call arguments were: $@"
  echo "You can rerun the delta/summarize.sh without reruning the model run(s). How-to: run ./benchmarks without arguments"
  echo "to get a help message that explains how to run the delta/summarize.sh without reruning the model."
  exit 1
fi

echo "Summary-processing: Deltareports... Arguments = $@"
echo "$0"
echo "The summarize script's pwd -P = $(pwd -P)"
cd ..
"$(dirname "$0")/${toolscript}" -c $@ | tee -a delta/$resultfile
