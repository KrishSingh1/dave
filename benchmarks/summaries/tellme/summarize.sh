#!/bin/bash

# This script is called with 3 arguments:
#   $1:	Name-tag for the new output
#   $2:	Name-tag for the comparison reference benchmark
#   $3:	Path to the comparison reference benchmark output
# The call to the Rscript tellme_global.Rmd here should have the folllowing syntax (email from Matt 2023-03-13)
# (3 of the 6 arguments are generated within this script).
#   tellme_global.Rmd <path to new benchmarking runs> <new name> <path to old run benchmarking runs> <old name> 
#   <path to evaluation data> <path to land use data>

# These paths are compute center specific. Update them to fit the center where you run your model.
# TODO - at LUnd need to make a shortcut on simba in /data/benchmark_data
# ln -s /data/evaluation_data/ Tellus
EVALDATA_PATH=${LPJG_BENCHMARK_DATA}/Tellus
LUDATA_PATH=${LPJG_BENCHMARK_DATA}/2023_03_02/landuse/LUH2/lu_1901_2015_luh2_Hist_CMIP_UofMD_landState_2_1_h_halfdeg_nourban_2019_11_15.txt


DEBUG_LOG_NAME="summarize.debug.log"

NAMETAG_NEW="$1"
NAMETAG_REF="$2"
NEW_OUTPUT_PATH="$(dirname $(pwd))"
REF_OUTPUT_PATH="$3"

# make some arguments based on location in filesystem and username
FULL_USERNAME=$(pinky -lb $(whoami) | cut -d: -f3 | tr -s " " | head -1 | sed -e 's/^[[:space:]]*//')
YAMLSETTINGS="tellme_global_config.yml"			# Do not change, unless you change also further down in this script where tellme_global.Rmd is hardcoded.
YAMLSETTINGS_PATHFILE="$(dirname $0)/${YAMLSETTINGS}"

if [ $# -ne 3 ]; then
  echo "Error: the summaries/tellme/summarize.sh script must be called with exctly 3 arguments." | tee $DEBUG_LOG_NAME
  echo "The call arguments were: $@"
  echo "You can rerun the tellme/summarize.sh without reruning the model run(s). How-to: run ./benchmarks without arguments"
  echo "to get a help message that explains how to run the tellme/summarize.sh without reruning the model."
  exit 1
fi

{
date | tee -a $DEBUG_LOG_NAME
echo "Tellme tool" | tee -a $DEBUG_LOG_NAME
echo "./benchmarks arguments = $@"   | tee -a $DEBUG_LOG_NAME
echo "Arguments to the tellus R-script:" | tee -a $DEBUG_LOG_NAME
echo "<path to new benchmarking runs> <new name> <path to old run benchmarking runs> <old name> <path to evaluation data> <path to land use data> <path to YAML file> <author>" | tee -a $DEBUG_LOG_NAME
echo "Argument 1 = $NEW_OUTPUT_PATH" | tee -a $DEBUG_LOG_NAME
echo "Argument 2 = $NAMETAG_NEW"     | tee -a $DEBUG_LOG_NAME
echo "Argument 3 = $REF_OUTPUT_PATH" | tee -a $DEBUG_LOG_NAME
echo "Argument 4 = $NAMETAG_REF"     | tee -a $DEBUG_LOG_NAME
echo "Argument 5 = $EVALDATA_PATH"       | tee -a $DEBUG_LOG_NAME
echo "Argument 6 = $LUDATA_PATH"     | tee -a $DEBUG_LOG_NAME
echo "Argument 7 = $YAMLSETTINGS_PATHFILE"     | tee -a $DEBUG_LOG_NAME
echo "Argument 8 = $FULL_USERNAME"     | tee -a $DEBUG_LOG_NAME
echo "'$0' = $0"      | tee -a $DEBUG_LOG_NAME
echo "pwd = $(pwd)"   | tee -a $DEBUG_LOG_NAME
echo | tee -a $DEBUG_LOG_NAME
} >/dev/null	# Don't output debug info to the console / slurm.log.


# Load the required software
if [[ "$ARCH" == "lunarc" ]]; then
module purge &>/dev/null
module load foss/2022a netCDF/4.9.0 CMake/3.23.1 Ghostscript/9.56.1
module load GCC/11.3.0  OpenMPI/4.1.4  R/4.2.1  Pandoc/3.1.2
elif (( "$ARCH" == "levante" )); then
    module load r/4.1.2-gcc-11.2.0 ghostscript/9.54.0-gcc-11.2.0
fi    

set -x		# Debug. Remove later.

# Here the call to the R-script tellme.Rmd
# Before call, the script is renamed to reflect what is being compared, and the users full name (Gecos) is inserted.
RSCRIPT="tellme_global.Rmd"			# Do not change, unless you change also further down in this script where tellme_global.Rmd is hardcoded.
RSCRIPT_PATHFILE="$(dirname $0)/${RSCRIPT}"
RSCRIPT_HTML="tellme_global.${NAMETAG_NEW}.vs.${NAMETAG_REF}.Rmd"
# MF: I would prefer not to do this cat/sed command, see instead the last argument below
cat $RSCRIPT_PATHFILE | sed "s/%%USER%%/$FULL_USERNAME/" >./$RSCRIPT_HTML
Rscript -e "rmarkdown::render('${RSCRIPT_HTML}',params=list(sim_directory=\"$NEW_OUTPUT_PATH\",sim_name=\"$NAMETAG_NEW\",reference_sim_directory=\"$REF_OUTPUT_PATH\",reference_sim_name=\"$NAMETAG_REF\",data_directory=\"$EVALDATA_PATH\",land_cover_file=\"$LUDATA_PATH\",config_file=\"$YAMLSETTINGS_PATHFILE\",author=\"$FULL_USERNAME\"))"

