#!/bin/bash

# Guess_diff. Diff files from LPJ-GUESS benchmarks as compared with a chosen reference run of benchmarks.
# Run script in the directory containing the benchmarks output folders (see variable bms below)
# Run as e.g. diff_outfiles_tslices.sh 5508 '/scratch/johan/Benchmarks/release_4.0/trunk_5508/output5508_all_pure'
# i.e. 	$1: a label naming the reference
#     	$2: path to the reference benchmarks output folders
# Output: is placed in each bm folder
# Metadata: this script is copied to pwd.
#	    $reflabel $refdir added to a file on path $0 for easy recollection of previous reference data
# Written by Johan Nord, 2019.


# Help message (when called with no arguments): shows top of this script
if [ $# -eq 0 ]; then
  echo
  echo Help text for $0
  tail -n+3 $0 | cut -c2- | head -9
  echo
  exit 0
fi

# This differs from the stand-alone version of this script.
# Retrieve the outputfolder, and then remove it from the argument list, so that the argument list is the same as in the stand-alone verion of this script.
outputfolder="$1"
shift

bms='crop_mixed_sites diurnal_pristine_sites emdi_europe emdi_global fluxnet europe pristine_sites secondary_stands global crop_global panarctic soil_temperature wetland_sites wetland_global tellus'

if [ $1 == "-1" ]; then
  bms=$2		# Only one benchmark, and which is specified by the user on the commandline
  reflabel=$3
  refdir=$4
else
  reflabel=$1				# E.g trunk8877
  [ -d $2 ]  || { ls -d $2; exit 1; }	# Test if reference dir exists, and exit if not.
  refdir=`cd $2; pwd`			# Change 2020-12-10	# E.g. ../output8877full, or '/scratch/johan/Benchmarks/release_4.0/trunk_5508/output5508_all_pure'
fi


# Define report file names
reportfile_out=diff_outfiles_vs_$reflabel.cout
reportfile_tslice61=diff_tslices1961to1990_vs_$reflabel.cout
reportfile_tslice91=diff_tslices1990to2000_vs_$reflabel.cout


# Error messages

if [ -e $reflabel ]; then
  echo 1st argument should be reflabel which must not match filename, 2nd arg path to reference benchmarks
  exit 99
fi
if [ ! -e $refdir ]; then
  echo 1st argument should be reflabel which must not match filename, 2nd arg path to reference benchmarks
  exit 99
fi

for bm in $bms; do	# Prevent overwriting of existing report files

 if [ -e $bm ];
 then
  if [ -e $refdir/$bm ];
  then

    cd $bm  >/dev/null

    if [ -f $reportfile_out ]; then
      echo $bm/$reportfile_out "already exists. Will be overwritten."
    fi

    if [ -f $reportfile_tslice61 ]; then
      echo $bm/$reportfile_tslice61 "already exists.  Will be overwritten."
    fi

    if [ -f $reportfile_tslice91 ]; then
      echo $bm/$reportfile_tslice91 "already exists.  Will be overwritten."
    fi

    cd - >/dev/null

  fi
 fi
done


# Main code

pwd -P

# Do the diffs

for bm in $bms; do

 if [ -e $bm ];
 then

  if [ -e $refdir/$bm ];		# Added 180627
  then

   cd $bm >/dev/null
   echo $bm
   #echo $refdir/$bm

   # tslices
   ls *1990.txt &> /dev/null && for file in $(ls *1990.txt); do diff -qs $file $refdir/$bm/$file; done > $reportfile_tslice61
   ls *2000.txt &> /dev/null && for file in $(ls *2000.txt); do diff -qs $file $refdir/$bm/$file; done > $reportfile_tslice91

   # .out files
   for file in $(ls *.out); do diff -qs $file $refdir/$bm/$file; done > $reportfile_out

   cd - >/dev/null

  else
   echo "Reference not found. Skipping\: $bm"        # Added 180627
  fi

 else
   echo "Not found. Skipping\: $bm"
 fi

done


# Report non-identical files

echo "Result summary"

{
echo "Diff against: $1 $2"
echo "COUNT NUMBER OF NON-IDENTICAL FILES PER EACH BM i.e. total files (T), identical files (I), files NOT identical (X)".

{
echo "${reportfile_out} T I X"
for f in $(ls */${reportfile_out}); do
  echo -n "$(dirname "$f")"
  echo -n " $(wc -l $f | cut -d\  -f1)"
  echo -n " $(grep -c identical $f)"
  echo -n " $(grep -c -v identical $f)"
  echo
done
} | column -t

{
echo "${reportfile_tslice61} T I X"
for f in $(ls */${reportfile_tslice61}); do
  echo -n "$(dirname "$f")"
  echo -n " $(wc -l $f | cut -d\  -f1)"
  echo -n " $(grep -c identical $f)"
  echo -n " $(grep -c -v identical $f)"
  echo
done
} | column -t

{
echo "${reportfile_tslice91} T I X"
for f in $(ls */${reportfile_tslice91}); do
  echo -n "$(dirname "$f")"
  echo -n " $(wc -l $f | cut -d\  -f1)"
  echo -n " $(grep -c identical $f)"
  echo -n " $(grep -c -v identical $f)"
  echo
done
} | column -t

} | tee "${outputfolder}/diff.${reflabel}.summary"

echo "LIST NON-IDENTICAL .OUT FILES. Only if it is not empty."
if [ -n "$(grep -v identical */$reportfile_out)" ]; then
  grep -v identical */$reportfile_out | sed "s/and//g" | sed "s/differ//g" \
  | sed "s/Files//g" | sed "s,/$reportfile_out,,g" | sed "s,: ,/,g" | tee "${outputfolder}/diff.${reflabel}.not-identical.list"
fi
