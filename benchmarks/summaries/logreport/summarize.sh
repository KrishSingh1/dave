#!/bin/bash

# Call guesslog_errors_summary.sh
# Writes a summary of the guess.logs
# Standard output is removed, to show only unexpected output.
# By Johan Nord, 2025.

toolscript="guesslog_errors_summary.sh"

echo "Summary-processing: generating log report... Arguments (should be empty) = $@"
echo "$0"
echo "The summarize script's pwd -P = $(pwd -P)"

if [ $# -ne 0 ]; then
  echo "Warning: the summaries/logreport/summarize.sh script doent use anu arguments."
  echo "The call arguments were: $@"
  echo "These arguments were ignored."
  echo
fi

"$(dirname "$0")/${toolscript}" ..	# Call it to run on output folders in one level up fron the folder "logreport".
