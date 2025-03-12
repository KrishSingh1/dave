#!/bin/bash

# Diff files from LPJ-GUESS benchmarks as compared with a chosen reference run of benchmarks.
# Run script in the directory containing the benchmarks output folders (see variable bms below).
# Example:
#   guess-diff-outfiles-tslices 5508 '/scratch/johan/Benchmarks/release_4.0/trunk_5508/output5508_all_pure'
# Argumnets:
#   $1: a label naming the reference
#   $2: path to the reference benchmarks output folders
# Output: is placed in each bm folder.
# Metadata: this script is copied to pwd.    !!!!!


toolscript="guess_diff.sh"
toolfoldername="$(basename "$(dirname "$0")")"
resultfile="${toolfoldername}/guess-diff-vs-${1}.summary.log"

if [ $# -ne 2 ]; then
  echo "Error: the summaries/diff/summarize.sh script must be called with exctly 2 arguments."
  echo "The call arguments were: $@"
  echo "You can rerun the diff/summarize.sh without reruning the model run(s). How-to: run ./benchmarks without arguments"
  echo "to get a help message that explains how to run the diff/summarize.sh without reruning the model."
  exit 1
fi

echo "Summary-processing: Output diff: are the outputs identical? Arguments = $@"
echo "$0"
echo "The summarize script's pwd -P = $(pwd -P)"
cd ..

# This call to guess_diff.sh" has one more argument that the stand alone version: the summary output-folder  name.
"$(dirname "$0")/${toolscript}" ${toolfoldername} $@ | tee -a $resultfile
