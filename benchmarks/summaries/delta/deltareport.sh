#!/bin/bash

# Comparison of data in the benchmark report catalog with comparator benchmark.
# By Johan Nord, 2017.
# Run in the BM outputfolder where there is a report directory to be compared.
# Output to new catalog report_delta.
# Input parameters:
# $1: suffix for the report_delta folder: e.g. trunk8538 => report_delta_trunk8538/
# $2: path to where the BM catalog is sitting,  where the report to compare with resides.

# Remains
# - Check att column headers är identiska

set +e

module purge &>/dev/null
module load GCCcore/11.3.0 Python/2.7.18 ImageMagick/7.1.0-37

summarytool="delta"
sub_imgs_cmd="$(dirname "$0")/deltareport_sub_imgs.sh"

# Help message (when called with no arguments): shows top of this script
if [ $# -eq 0 ]; then
  echo
  echo Help text for $0
  tail -n+3 $0 | head -n5       # Change head -n to the number of helptext lines
  echo
  exit 0
fi


### Definitions

suffix=$1
path_comparisonrun=$2

replacementpath_to_catalog_bmcommon="deltareport_altsource_bmcommon"


### functions

padfloat() {

   for x in "$@"
   do
     { printf "%14.2f " $x || true ; }
   done
}

padstring() {

   for x in "$@"
   do
#echo padstring $X
      printf "%14s " $x
   done
}

# Test if string is a perfect numeric
testnumeric() {

  string=$1   # GÖR LOCAL!!!

  # Trim EXTERNAL_SPACE           # källa: https://stackoverflow.com/questions/369758/how-to-trim-whitespace-from-a-bash-variable
  string="$(echo -e "${string}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
  # debug: echo $string

  # Allow leading negative
  string="$(echo -e "${string}" | sed -e 's/^\-//')"

  # Allow decimalperiod - !!!: Weak test - does not care about position: only within string...
  string="$(echo -e "${string}" | sed -e 's/\.//')"

  # Allow numerics		# Källa: https://stackoverflow.com/questions/806906/how-do-i-test-if-a-variable-is-a-number-in-bash
  case $string in
    ''|*[!0-9]*) echo "" ;;	# return false
    *) echo $string ;;		# return true
  esac

}

# Return n:th line
# $1 = line number
# $2 = file
tailhead1() {

  tail -n+${1} ${2} | head -n1

}

### Main code ###


### Init

#set -e

reportpath=$(pwd -P)/report
diffreportpath=$(pwd -P)/report_delta41_$suffix
if [ ! -d $reportpath ]; then
  echo "No such path: " $reportpath
  exit 1
fi

benchmark=$(basename $(pwd -P))
comparatorpath=`cd $2/${benchmark}/report; pwd -P`
if [ ! -d "$2/${benchmark}/report" ]; then
  echo "No such path: " "$2/${benchmark}/report"
  exit 1
fi

# Diagnostic

pwd -P
ls -d $reportpath
echo Generating delta report based on following comparator benchmark
ls -d $comparatorpath
mkdir -p $diffreportpath   #-p: no error if exists
if [ ! -d $diffreportpath ]; then
  echo "No such path: " $diffreportpath
  exit 1
fi

cd $diffreportpath
ls -d $diffreportpath

# Make a reference to the comparator src
rm -f comparator_benchmark
ln -s $comparatorpath comparator_benchmark

# Make a reference to the target of symlink common for purpose of reading script code
# and if the original svn code is not present, make a symlink to a replacement svn code.
# Later, improve the script by, if missing, having it svn-download it using info.txt info.
echo "PWDP"			# DEBUG
pwd -P				# DEBUG
if ls -q ../common/; then
  commonlinkpath=`cd ..; cd common; pwd -P`
  rm -f srccatalog_bmcommon
  ln -s $commonlinkpath srccatalog_bmcommon
else
  echo The symlink common is broken...
  if [ -e "../../${replacementpath_to_catalog_bmcommon}/." ]; then
    echo "Making symlink to replacement code as per ${replacementpath_to_catalog_bmcommon}"
    commonlinkpath=`cd ../..; cd ${replacementpath_to_catalog_bmcommon}; pwd -P`
#    ln -s `cd ../..; cd ${replacementpath_to_catalog_bmcommon}; pwd -P` srccatalog_bmcommon
    rm -f srccatalog_bmcommon
    ln -s $commonlinkpath srccatalog_bmcommon
  else
    echo "Alt common: No such symlink: ../../${replacementpath_to_catalog_bmcommon}. Aborting."
    echo "To proceed, make a symlink on the level of info.txt called ${replacementpath_to_catalog_bmcommon}"
    echo "which points to the benchmarks/common catalog if the relevant LPJ-GUESS source code."
    exit 1
  fi
fi

#set +e   ## Debug 190308


###################################################
### Convert tables into comparison tables


# Loop though report data files

infiles=$(grep 'textfile src' ${reportpath}/report.xml | cut -f2 -d\")

for infile in $infiles; do

  # If output as txt instead of html is required
  echo >> diffreport.txt
  echo $infile >> diffreport.txt

  #
  filerep=${reportpath}/$infile
  filecomp=${comparatorpath}/$infile
  filediff=${diffreportpath}/$infile

  #################
  # Determine what kind of table this is: 1) 1-data-row table (aslicetype), 2) vertical table, 3) multirow normal table (horizontal).
  # The type-2 test overrides the two other types.

  tabletype=1	# Default type if not found below to be a differnet type (2 or 3).

  # The table is deemed vertical (2) if in second row any value is nonnumeric,
  # and if not that, it is deemed multirow (3)  if it has more than 2 lines and is NOT vertical type (2).

  # Test for type multirow (3)
  tablenumberoflines=$(wc -l $filerep | cut -f1 -d\  )
  two=2
  #echo $filerep tablenumberoflines $tablenumberoflines
  if [ ${tablenumberoflines} -gt $two ]; then 
    if [ ${tablenumberoflines} -lt $two ]; then 
      echo Wrong table? Empty file? Only one line - is that ok?
    fi
    tabletype=3
    #echo Multirow table 	# debug
  fi

  # Test for type vertical (2)
  # Loop throu second row (i.e. the line after the header line) and test each item for numericallity.
  secondrow=$(cat $filerep | tail -n+2 | head -n1)

  for col in $secondrow; do
    if [ -z $(testnumeric $col) ]; then
      tabletype=2
      #echo Vertical table. col\= $col	 # debug
    fi
  done
  #
  #
  ################

  # Debug information. Table type
  echo Tabletype conclusion: $infile Tabletype \= $tabletype
  #  1 Std single row table
  #  2 Vertical table
  #  3 Multirow table

  ## If type-2 table (vertical):   ###################
  if [ ${tabletype} -eq 2 ]; then
set +e
    paste ${filerep} ${filecomp} >$filediff

    cat $filediff >> diffreport.txt
#set -e

  else
  ## If type-1 and -3 table:   ###################
  ## Compute differences and make comparison-lines items


  # NOT USED Count number columns in report file
  ncols=$(head -1 ${reportpath}/$infile | wc -w)

  # Get list of columnheaders from report file
  cols=$(head -1 $filerep)

  # Rename columnheaders in comparator file, but keep them intact in list compcols
  compcolheader_suffix=2		#_cmpxyzzyx
  compcols=$(head -1 $filecomp)
  {
    for col in $compcols; do
      echo -n "${col}${compcolheader_suffix} "
    done
    echo
  } > filecomp_rename.tmp
  tail -n+2 $filecomp >>filecomp_rename.tmp

  # Copy table data to working file (report file, and header-modified comparator file)
  # Use only lines down to before an empty line. Lines after an empty line are considered to be commenting.
  paste ${filerep} filecomp_rename.tmp | awk 'NF { print } !NF { exit }' >pastefile.tmp


  # Compute and make comparison-lines items

  echo -n "compute pastefile.tmp -n -i " >compute_absdiff.sh.tmp
  echo -n "compute pastefile.tmp -n -i " >compute_reldiff.sh.tmp

  echo -n >nc.tmp
  for i in $(seq -s ' ' 2 ${tablenumberoflines}); do
    echo "n.c." >> nc.tmp
  done

  # Operation: a whole column at a time; same for both type-1 and type-3 tables.
  echo -n >filecomp.tmp
  different_columnheaders=""
  i=1
  for col in $cols; do
     if [[ $compcols == *"$col"* ]] ; then		# i.e. in pseudocode: if [ col in $compcols ]; then

	compcol=${col}${compcolheader_suffix}
	echo -n "'${col}-${compcol}' " >>compute_absdiff.sh.tmp
	echo -n "'(${col}-${compcol})/${compcol}*100' " >>compute_reldiff.sh.tmp

	cat $filecomp | tr -s " " | cut -d\  -f2- | cut -d\  -f$i | paste filecomp.tmp - > filecomp.tmp2
	mv filecomp.tmp{2,}
	echo $i $col $compcol        # debug

     else        # i.e. if column header is not existing in comparator table then:

        different_columnheaders="1"
	echo -n "'0' " >>compute_absdiff.sh.tmp
	echo -n "'0' " >>compute_reldiff.sh.tmp

        echo "$col" >compnc.tmp
        cat nc.tmp >>compnc.tmp
	cat $filecomp | tr -s " " | cut -d\  -f2- | cut -d\  -f$i | paste compnc.tmp - > filecomp.tmp2
	mv filecomp.tmp{2,}
	echo $i $col nc $compcol	# debug

     fi
     ((i+=1))
  done
  cp filecomp{,.final}.tmp

  # Do the calculation
  echo -n "-o" absdiff.tmp >>compute_absdiff.sh.tmp
  echo -n "-o" reldiff.tmp >>compute_reldiff.sh.tmp
  echo  >>compute_absdiff.sh.tmp
  echo  >>compute_reldiff.sh.tmp
  sh ./compute_absdiff.sh.tmp &>compute.log
  sh ./compute_reldiff.sh.tmp &>compute.log


  ## Make the data files for the new xml-file

{
  padstring $(head -n1 $filerep); echo
  if [ ! -z $different_columnheaders ]; then
    padstring $(head -n1 $filecomp); echo
  fi

  for i in $(seq -s ' ' 2 ${tablenumberoflines}); do

    padfloat $(tailhead1 $i $filerep			); echo "          "Report
    padfloat $(tailhead1 $i filecomp.final.tmp 	); echo "          $suffix"
    padfloat $(tailhead1 $i absdiff.tmp 		); echo "          "Difference
    padfloat $(tailhead1 $i reldiff.tmp 		); echo "          "Rel diff %;

    if [ ! "$i" -eq "${tablenumberoflines}" ]; then
      echo " "
    fi

  done

} >$filediff

  cat $filediff >> diffreport.txt

fi	# type of table = 1 and 3

done	# next table file


# Clean up tmp files
#rm compute_absdiff.sh.tmp compute_reldiff.sh.tmp
#rm absdiff.tmp reldiff.tmp pastefile.tmp compute.log


################################
### Make the new html file


# Make the new xml-file
cp $reportpath/report.xml $diffreportpath

# Make the combined svn-info file
{ cat $reportpath/../../info.txt; echo; echo Comparator benchmark; echo; echo; cat $comparatorpath/../../info.txt; } >svn.info.txt

# run report2html
python2 `cd srccatalog_bmcommon; pwd -P`/../report2html report.xml svn.info.txt >index.html

# Compar-info in header
sed -i 's|</p><table|\n</p><table|' index.html
{
  echo -n $(head -n1 ${diffreportpath}/index.html)	   # Do all inserts after what is here made to be line one in html script (not in what-you-see)

  echo -n ' <p\><h2>Comparison delta report:</h2></p>'

  echo -n Present benchmark\ 
  echo -n $(grep "VCS Path" ${diffreportpath}/index.html | tail -n+1 | head -n1 | cut -d\  -f3)
  echo -n \ '<a href="../report/index.html">'
  echo -n r$(grep "Revision" ${diffreportpath}/index.html | tail -n+1 | head -n1 | cut -c11-)
  echo -n '</a>'

  echo -n \ versus reference\ 
  echo -n $(grep "VCS Path" ${diffreportpath}/index.html | tail -n+2 | head -n1 | cut -d\  -f3)
  echo -n \ "<a href="file://$comparatorpath/index.html">"
  echo -n r$(grep "Revision" ${diffreportpath}/index.html | tail -n+2 | head -n1 | cut -c11- )
  echo -n '</a>'" ($suffix)."

  echo -n '<br>'
  echo -n 'Tables: row 1: present run, row 2: reference, row 3: difference, row 4: relative difference (%).'
  echo -n '<br>'
  echo -n 'Map graphs: Present run is to the left, reference benchmark in the middle (or to the right), difference plot to the furthest right (no difference plot  means the two other plots are identical)'
  echo -n '<br>'

  tail -n+2 ${diffreportpath}/index.html		# And finally add all of the rest if the file
} > ${diffreportpath}/index.tmp
mv ${diffreportpath}/index.tmp ${diffreportpath}/index.html



### Finish

cd ..
echo
echo Finished tabelling $(basename $0) $2
echo now doing images and reformatting html code ...

bash $sub_imgs_cmd $diffreportpath	#N.B. is a path, to only a folder name     # upvote! https://stackoverflow.com/questions/6121091/get-file-directory-path-from-file-path

echo Finished images and reformatting html code. All finished.
