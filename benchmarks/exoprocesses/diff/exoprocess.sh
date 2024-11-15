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


resultfile="guess-diff.vs-${1}.result"

echo "Exo-processing: Output diff: are the output identical... Arguments = $@"
echo "$0"
echo "The exoscript's pwd -P = $(pwd -P)"
cd ..
guess-diff-outfiles-tslices $@ | tee -a $resultfile
