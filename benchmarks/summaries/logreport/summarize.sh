#!/bin/bash

# Summmary tool logreport. It writes a summary of the guess.logs of all the benchmarks runs.
# Standard output is removed, to show only unexpected output.
# This script calls the main script guesslog_errors_summary.sh
# By Johan Nord, 2022.
# Usage
# ./benchmarks -i "fluxnet tellus" -s "logreport [-h] <path to reference output>" <benchmark-runs output directory>
# Using option -h will cause logreport to print this help message and exit.

toolscript="guesslog_summary.sh"

if [ $# -ne 1 ]; then
  echo "Error: the summaries/logreport/summarize.sh script needs one argument:"
  echo "The path to the directory with reference output, e.g. from latest trunk."
  echo "The call arguments were: $@"
  echo "Abort."
  echo
  HELP="1"		# Help will cause exit 1
fi

# Handle the command line arguments
while getopts ":h" opt; do
    case $opt in
        h ) HELP="1" ;;
		\? ) echo "Error! Unavailble option $opt. Abort!"; HELP="1" ;;		# Help will cause exit 1
    esac
done
shift $((OPTIND-1))

# Display help msg and exit
if [ "$HELP" == "1" ]; then
  echo Help text for $0
  echo
  tail -n+3 $0 | cut -c2- | head -n7        # Change head -n to the number of helptext lines
  echo
  exit 1
fi

# Main code

echo "Summary-processing: generating log report..." 
echo "Arguments (should be 1, the path to reference log files) = $@"
echo "$0"
echo "The summarize script's pwd -P = $(pwd -P)"

# Call the main script
"$(dirname "$0")/${toolscript}" $1 ..	# Call it to run on output folders in one level up fron the folder "logreport".
