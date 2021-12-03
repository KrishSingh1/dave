#!/bin/bash

echo "Tellme tool" | tee postpostprocess.dummy-output
echo "Arguments = $@" | tee -a postpostprocess.dummy-output
echo "Argument 1 = $1" | tee -a postpostprocess.dummy-output
echo "Argument 2 = $2" | tee -a postpostprocess.dummy-output
echo "pwd -P = $(pwd -P)" | tee -a postpostprocess.dummy-output

