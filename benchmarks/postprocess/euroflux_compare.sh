#!/bin/bash

# This script produces a scatter plot for comparison of two columns from
# an input file. It is meant to be used with the output from euroflux
# benchmarks, so the following applies:
#
# 1. If a value in one of the chosen columns is -9999.0, the line is ignored.
# 2. The first line is used for titles for the x and y axes

if [ $# -lt 5 ]; then
    echo "Usage: $0 <filename> <column1> <column2> <filename>.png <title>"
    echo
    echo "Columns are 1-based (1 = first column)."
    echo "column1 is used for the y-axis, column2 for the x-axis."
    echo "For instance:"
    echo "$0 eurofluxmonthly.out 6 5 euroflux_nee.png \"NEE\""
    exit 1
fi

source `dirname $0`"/scatter_plot.sh"

# Get the arguments
DATA_FILE=$1
COLUMN1=$2
COLUMN2=$3
OUT_FILE=$4
TITLE=$5

# Create a temporary file for only the data to plot
TEMP_FILTERED=$(mktemp)

# Extract the chosen columns for all lines where neither value is -9999.0,
# and skip the first line
CONDITION='NR>1 && $'$COLUMN1' != -9999.0 && $'$COLUMN2' != -9999.0'
ACTION='{print $'$COLUMN1'" " $'$COLUMN2'}'
awk "${CONDITION}"' '"${ACTION}" ${DATA_FILE} > ${TEMP_FILTERED}

# Extract the xlabel and ylabel
XLABEL=$(awk 'NR==1 {print $'${COLUMN1}'}' ${DATA_FILE})
YLABEL=$(awk 'NR==1 {print $'${COLUMN2}'}' ${DATA_FILE})

# Call the scatter_plot function to do the actual plotting
scatter_plot "${TITLE}" "${XLABEL}" "${YLABEL}" ${TEMP_FILTERED} ${OUT_FILE}

# Clean up temp files
rm ${TEMP_FILTERED}
