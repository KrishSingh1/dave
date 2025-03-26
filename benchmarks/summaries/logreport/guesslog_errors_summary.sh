#!/bin/bash
# guesslog_errors_summary (standallone version is guess-grep-guesslog-errors-allBMs.sh aka guess-log-error-extract)
# By Johan Nord, 2025.
# Output is written to a file, see variable outfile below. It is not printed to std out.
# Arguments:
# If 1 argument, operates on current dir. Else (i.e. 2 arguments) it operates on dir given in $2.
# $1: dir containing bm folders with guess.log files.
# $2: (optional) path to operate on. Only used when used stand-alone, not togeher with logreport benchmarks summary.


### Special tweaks

# Exclude BM(s?)
EXCLUDEDBMS=""
#"fluxnet"
#"tellus"

# Exclude these words or strings, separate with |
FILTERPATTERN="C pool change|C flux|Period C balance|C balance year|N pool change|N flux|Period N balance|N balance year"	# elapsed|^Commen|^Last|Using|C pool change|C fluxes|Period C|^$
# and do it only on these BMs:
FILTERBMS="tellus"

bms="crop_global crop_mixed_sites diurnal_pristine_sites emdi_europe emdi_global europe fluxnet global panarctic\
 pristine_sites secondary_stands soil_temperature tellus wetland_global wetland_sites"


### Initialisation

# This is where output is written to. It is not printed to std out.
# Then cd to where the log files analusis is done
outfile="`pwd -P;`/guesslog_errors.txt"
if [ $# -ne 0 ]; then
  cd $2			
fi

# Path to the reference output 
expected_files_reference="$1"
if [ ! -d "$expected_files_reference" ]; then
   echo "The path to the reference output is invalid. Abort."
   exit 1
fi


### Main code


### Intro
string_intro=$(
#echo $0     #> $outfile
echo "Analysis of:"
pwd -P      #>> $outfile
echo
echo "This output is also saved in in file $outfile."
echo
echo "1a. Look for missing report/index.html files."
echo "1b. Look for empty .out files."
echo "2.  Print only unexpected output."
echo "3.  Counting start- and finished-lines. The numbers should match."
echo "4.  Check that number of outputfiles are ok (not missing), i.e.:"
echo "    checking the presence of expected files or superfluous files as cmp w a BM reference."
echo "5.  All output except expected warnings, i.e. start- and finished-lines, plus unexpected output."
)


### Verify that report-folders have been created
string_missingreport=$(
echo
echo "_________________________________________________________"
echo "Look for missing report/index.html. This should be empty."
diff -u <(cd ${expected_files_reference};\
 ls -d1 */report/index.html) <(ls -d1 */report/index.html)\
 | grep -v ^\-\-\- | grep -v ^+++ | grep "^\-\|^+" | less
echo
)
#} > $outfile 2>&1


### Look for empty outfiles
string_emptyoutfiles=$(
echo
echo "_______________________________________________________________________________________________"
echo "Look for empty outfiles. This should be empty. This is not not a comparison with the reference."
wc -l */*.out | awk '$1==0'
echo
)


### Print guesslog errors and finished-lines from logs
{	#string_guesslogerrors=$(
echo
for bmcat in $bms; do
  if [[ $EXCLUDEDBMS =~ $bmcat ]]; then
    echo "Excluded $bmcat"
    echo
    continue
  fi
  echo $bmcat
  GENERAL="^(LPJ-GUESS cohort mode - |Using soil code and Nitrogen deposition for|LPJ-GUESS test Secondary stands|Description: |Last year of cropland fraction data used from year 2016 and onwards)"
  if [[ $FILTERBMS =~ $bmcat ]]; then
    FILTER="${GENERAL}|${FILTERPATTERN}"
    echo "Extra excluded lines by these filter terms: $FILTERPATTERN"
  else
    FILTER="${GENERAL}"
  fi

  # Remove standard output
  grep -v 'Commencing simulation' $bmcat/guess.log \
  | grep -E -v "complete.*elapsed.*remaining" | grep -v -e '^[[:space:]]*$' \
  | grep -v 'fraction data used from year 2007 and onwards' | grep -v "LPJ-GUESS cohort mode \- global pfts" \
  | grep -E -v "\[LPJ-GUESS  .*2017\]" | grep -v "\-\-\-\-\-\-\-\-\-\-" \
  | awk '{\
    if ( $0 ~ /^~~~~~~~~~~~~/ ) {\
      if ( skipline==0 ) { skipline=1 } else { skipline=0 } } else { if ( skipline==0 ) { print $0 } }\
    }' \
  | grep -vE "$FILTER"		# note the backslash at the end of the prev line.
  echo

done

} > $outfile 2>&1
string_guesslogerrors=$(cat $outfile)


string_unexpected=$(
echo
echo "______________________"
echo "Unexpected output only"

cat $outfile | grep -v Finished |grep -v "^\[LPJ-GUESS" \
| grep -v ^"LPJ-GUESS cohort mode - " \
| grep -v "^Using soil code and Nitrogen deposition for" \
| grep -v "^LPJ-GUESS test Secondary stands" \
| grep -v "^Description: " \
| grep -v "^Last year of cropland fraction data used from year 2016 and onwards"
# note the backslash at the end onf the prev line.

echo "________________________________________________"
echo "Unexpected output plus Start- and Finished-lines"

nLPJG=$(cat $outfile | grep -c "\[LPJ-GUESS")
nFinished=$(cat $outfile | grep -c "Finished")

echo "Number of \[LPJ-GUESS found: $nLPJG"
echo "Number of Finished found: $nFinished"
echo
echo "______________________________________________"
echo "files as cmp w $expected_files_reference"
echo "(diff of ls of each bm-catalog found in . with the reference):"

for bmcat in $bms; do

  if [ -d $bmcat ]; then
    echo
    echo $bmcat
    diff -u <(ls "${expected_files_reference}/$bmcat") <(ls $bmcat) |\
         grep -e "^\-\|^\+"  |\
         grep -v slurm | grep -v diff_out | grep -v diff_t |\
         grep -v "report_delta" | grep -v '\-\-\- /dev/fd/' | grep -v '\+\+\+ /dev/fd/'
  fi
done

echo
)
#} >> $outfile 2>&1

echo "$string_intro" > $outfile 2>&1
echo >> $outfile 2>&1
echo "$string_missingreport" >> $outfile 2>&1
echo >> $outfile 2>&1
echo "$string_emptyoutfiles" >> $outfile 2>&1
echo >> $outfile 2>&1
echo "$string_unexpected" >> $outfile 2>&1
echo >> $outfile 2>&1
echo "__________________________________________________" >> $outfile 2>&1
echo "Print finished-lines from logs, and errors if any." >> $outfile 2>&1
echo "$string_guesslogerrors" >> $outfile 2>&1

echo
echo Finshed
echo The output is saved in $outfile
