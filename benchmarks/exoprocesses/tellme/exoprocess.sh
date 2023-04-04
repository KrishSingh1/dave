#!/bin/bash

# This script is called with 3 arguments:
# $1:	Name-tag for the new output
# $2:	Name-tag for the comparison reference benchmark
# $3:	Path to the comparison reference benchmark output

# The call to the Rscript tellme.Rmd here should have the folllowing syntax (email from Matt 2023-03-13)
# (3 of the 6 arguments are generated within this script).
# tellme.Rmd <path-to-new-benchmarking-runs> <new-name> <path-to-old-run-benchmarking-runs> <old-name> 
# <path-to-benchrmarking-data-dir> <path-to-land-use-data>

# These paths are compute center specific. Update to fit the center where you run your model.
DATA_PATH=/data/benchmark_data/
LUDATA_PATH=/data/benchmark_data/2019-04-29/landuse/


RSCRIPT="tellme.Rmd"
DEBUG_LOG_NAME="exoprocess.debug.log"

NAMETAG_NEW="$1"
NAMETAG_REF="$2"
NEW_OUTPUT_PATH="$(dirname $(pwd))"
REF_OUTPUT_PATH="$3"
RSCRIPT_PATHFILE="$(dirname $0)/${RSCRIPT}"


if [ $# -ne 3 ]; then
  echo "Error: the exoprocesses/tellme/exoprocess.sh script must be called with exctly 3 arguments." | tee $DEBUG_LOG_NAME
  echo "The call arguments were: $@"
  exit 1
fi


echo "Tellme tool" | tee $DEBUG_LOG_NAME
echo "./benchmarks arguments = $@"   | tee -a $DEBUG_LOG_NAME
echo "Arguments to the tellus R-script:" | tee -a $DEBUG_LOG_NAME
echo "<path-to-new-benchmarking-runs> <new-name> <path-to-old-run-benchmarking-runs> <old-name> <path-to-benchrmarking-data-dir> <path-to-land-use-data>" | tee -a $DEBUG_LOG_NAME
echo "Argument 1 = $NEW_OUTPUT_PATH" | tee -a $DEBUG_LOG_NAME
echo "Argument 2 = $NAMETAG_NEW"     | tee -a $DEBUG_LOG_NAME
echo "Argument 3 = $REF_OUTPUT_PATH" | tee -a $DEBUG_LOG_NAME
echo "Argument 4 = $NAMETAG_REF"     | tee -a $DEBUG_LOG_NAME
echo "Argument 5 = $DATA_PATH"       | tee -a $DEBUG_LOG_NAME
echo "Argument 6 = $LUDATA_PATH"     | tee -a $DEBUG_LOG_NAME
echo "'$0' = $0"      | tee -a $DEBUG_LOG_NAME
echo "pwd = $(pwd)"   | tee -a $DEBUG_LOG_NAME

# Here the call to the Rscript tellme.Rmd
$RSCRIPT_PATHFILE $NEW_OUTPUT_PATH $NAMETAG_NEW $REF_OUTPUT_PATH $NAMETAG_REF $DATA_PATH $LUDATA_PATH
