#!/bin/bash

# See comment in benchmarks/secondary_stands/guess.ins about this test and other output of this benchmark.

cat cmass_sts.out | tail -n+2 | head -n50\
 | awk '$4 != $5 {print "\nTEST FAILED! Benchmark test secondary_stands failed: Natural and mixed columns differ at line", NR, "in cmass_sts.out","\nThis indicates a problem with secondary stands (e.g. a member variable was not inititated in constructor, etc.)"; exit;}'\
 >testresult.txt

cat testresult.txt >>guess.log


# Benchmark report

describe_benchmark "LPJ-GUESS - Secondary stands technical test"

if [ -z "$(cat testresult.txt)" ]; then		# Test passed.
  echo "Test passed." >testresult.txt
fi						# If test failed, use the pre-existing text in testresult.txt

describe_textfile testresult.txt

