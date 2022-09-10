#!/bin/bash

# See comment in benchmarks/secondary_stands/guess.ins about this test and other output of this benchmark.

cat cmass_sts.out | tail -n+2 | head -n50\
 | awk '$4 != $5 {print "\nTEST FAILED! Benchmark test secondary_stands failed: Natural and mixed columns differ at line", NR, "in cmass_sts.out","\nThis indicates a problem with secondary stands (e.g. a member variable was not inititated in constructor, etc.)"; exit;}'\
 >testresult.txt

cat testresult.txt >>guess.log

if [ "$(wc -l testresult.txt)" != "0" ]; then	# Test failed
#  describe_benchmark "LPJ-GUESS - Secondary stands technical test \n$(cat testresult.txt)"
  echo failed $(wc -l testresult.txt)
else			# Test passed.
#  describe_benchmark "LPJ-GUESS - Secondary stands technical test \n\nTest passed."
  echo ok $(wc -l testresult.txt)
fi


#This does not work, unfortunately. Syntax error from awk at '!'"
#awkmsg1="TEST FAILED! Natural and mixed columns differ at line";\
#awkmsg2="\nThis indicates a problem with secondary stands";\
#cat cmass_sts.out | tail -n+2 | head -n51\
# | awk -v awkmsg1=$awkmsg1 -v awkmsg2=$awkmsg2 '$4 != $5 {print $awkmsg1, NR, "in cmass_sts.out",$awkmsg1; exit;}'\
# | less
