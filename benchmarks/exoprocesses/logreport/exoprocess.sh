#!/bin/bash

# Writes a summary of the guess.logs
# Standard output is removed, to show only unexpected output.


resultfile="guess-log-report.txt"

echo "Exo-processing: generating log report... Arguments (should be empty) = $@"
guess-log-error-extract | tee -a $resultfile

