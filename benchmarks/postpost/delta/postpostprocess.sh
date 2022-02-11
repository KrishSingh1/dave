#!/bin/bash

# Decision by Johan and Stefan 2022-02-11 
# to rename postpostprocess-files and folders to exoprocess.
# because postpostprocess is impossible to pronounce.
# exo meaning "on the outside" etc.
# e.g.:
#      benchmarks\postpost\delta\postpostprocess.sh"
#  --> benchmarks\exoprocesses\delta\exoprocess.sh"
#					^					^

echo "Delta tool" | tee postpostprocess.dummy-output
echo "Arguments = $@" | tee -a postpostprocess.dummy-output
echo "Argument 1 = $1" | tee -a postpostprocess.dummy-output
echo "Argument 2 = $2" | tee -a postpostprocess.dummy-output
echo "pwd -P = $(pwd -P)" | tee -a postpostprocess.dummy-output
