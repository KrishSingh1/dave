#!/bin/bash

# deltareports runs deltareport on all benchmarks folders in the current directory.
# Written Johan Nord, 2017
# deltareport in turn compares two LPJ-GUESS benchmarks reports - the current one and a
# reference report, and outputs a delta report in each benchmark's directory.
# Usage:
# Deltareports [-options] <label> <path to comparison reference data>
# options:
# -h	This help text
# -c	Use option -c when this script is run on a cluster node. With this option, it
#	will not run the deltareport jobs with nohup, which doesnt work well on a cluster.

module load Python/2.7.18
module load ImageMagick/7.1.0-37

summarytool="delta"
deltareportcmd="$(dirname "$0")/deltareport.sh"

# Options - Handle the command line switches
while getopts "ch" opt; do
    case $opt in
        c ) CLUSTER="1" ;;  	# If the script is un in a cluster it will not run the deltareport jobs under nohup, which doesnt work well on a cluster.
        h ) HELP="1" ;;         # Display help msg and exit
    esac
done
shift $((OPTIND-1))


# Help message which shows the top of this script
if [ "$HELP" == "1" ]; then
  echo
  echo Help text for $0
  tail -n+3 $0 | cut -c2- | head -10       # Change head -n to the number of helptext lines
  echo
  exit 0
fi


if [ ! -d "$2" ] || [ ! $# -eq 2 ] ; then
  echo "Comparator path not found, or wrong number of arguments. Exit."
  echo "This script takes 2 arguments: a label, and the path to the folder containing the benchmarks to compare with."
  exit 1
fi

label=$1
comparator=$(realpath -e $2)

tmpfile_reportsnotfound="$HOME/.deltareports.notfound.tmp"
tmpfile_reportsfound="$HOME/.deltareports.found.tmp"

echo -n "" >${tmpfile_reportsnotfound}  2>/dev/null
echo -n "" >${tmpfile_reportsfound}  2>/dev/null

echo "Comparison with: $comparator"
echo -n "Starting one process per benchmark..."

#bms=$( ls )
bms="emdi_europe emdi_global pristine_sites diurnal_pristine_sites fluxnet europe soil_temperature wetland_sites wetland_global panarctic crop_global global"

for bm in  $bms; do

 if [ -d $bm ]; then

  if cd $bm ; then
    if ls -d report &>/dev/null; then
      if [ "$CLUSTER" != "1" ]; then
        nohup "$deltareportcmd" $label $comparator &>${summarytool}_deltareport_${label}.ceout &
      else
        echo "Doing $bm" &>${summarytool}_deltareport_${label}.ceout
        "$deltareportcmd" $label $comparator &>>${summarytool}_deltareport_${label}.ceout		# I.e. don't use nohup on a cluster node.
      fi
      if [ -d "$comparator/$bm/report" ]; then
        echo -n " $bm" >>${tmpfile_reportsfound}  2>/dev/null
      else
        echo -n " $bm*">>${tmpfile_reportsnotfound}  2>/dev/null
      fi
    else
      echo -n " [$bm: no local report]">>${tmpfile_reportsnotfound}  2>/dev/null
    fi
    cd ..
  fi

 else
   echo -n " $bm" >>${tmpfile_reportsnotfound}  2>/dev/null
 fi

done

cat ${tmpfile_reportsfound}  2>/dev/null
echo
echo -n "All bms deltareport processes started"
if [ ! -z "$(cat "${tmpfile_reportsnotfound}" 2>/dev/null)" ]; then
  echo " except, not found:$(cat ${tmpfile_reportsnotfound})"  2>/dev/null
else
  echo
fi
echo "Each deltareport can be found in each benchnmark output's folder, in a folder report_delta_$label"
