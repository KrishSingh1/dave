dnl       M4 template for the submit script
dnl
dnl       The submit script will look exactly like this file
dnl       except that these comments are removed and the words
dnl       BINARY and DATE will be replaced as appropriate.
dnl
dnl       The changecom command below makes sure keywords are
dnl       replaced inside bash comments as well.
changecom()dnl
#!/bin/bash
#
# submit.sh
#
# Portable bash script to run LPJ-GUESS version:
# BINARY
# as a parallel job using PBS on Milleotto.
#
# Created automatically on DATE
# 
# Usage: 
#
#   1. Copy script to the directory where you want output written.
#      This will be called the RUN DIRECTORY.
#   2. In an editor, set appropriate values for the variables NPROCESS,
#      INSFILE, GRIDLIST and OUTFILES (NB: no space after the = sign):

NPROCESS=15
WALLTIME=150:00:00
INSFILE=guess.ins
GRIDLIST=gridlist.txt
OUTFILES='cmass.out firert.out anpp.out lai.out cflux.out dens.out tot_runoff.out mgpp.out mra.out mrh.out cpool.out mnpp.out mlai.out mnee.out maet.out mpet.out mevap.out mintercep.out mrunoff.out mwcont_upper.out mwcont_lower.out'

#      Where:
#      NPROCESS = number of processes in parallel job
#      WALLTIME = maximum wall (real) time for job hh:mm:ss
#      INSFILE  = path to ins file from run directory
#      GRIDLIST = path to gridlist file from run directory
#      OUTFILES = list of LPJ-GUESS output files in single quotes,
#                 and separated by spaces (filenames only, including
#                 extension, no directory.)
#
#   3. Run the script using the command:
#        sh submit.sh
#      or:
#        sh submit.sh [-n <name>] [-s <file>]
#
#      Both arguments are optional and interpreted as:
#      name     = the name of the job (shown in PBS queue)
#      file     = filename of a file which can override the variables
#                 above
#
# Nothing to change past here
########################################################################

# Handle the command line arguments
while getopts ":n:s:" opt; do
    case $opt in
	n ) name=$OPTARG ;;
	s ) submit_vars_file=$OPTARG ;;
    esac
done

# Override the submit variables with the contents of a file, if given
if [ -n "$submit_vars_file" ]; then
    source $submit_vars_file
fi

# This function creates the gridlist files for each run by splitting
# the original gridlist file into approximately equal parts.
function split_gridlist {
    # Create empty gridlists first to make sure each run gets one
    for ((a=1; a <= NPROCESS ; a++)) 
    do
      echo > run$a/$GRIDLIST
    done

    # Figure out suitable number of lines per gridlist, get the number of
    # lines in original gridlist file, divide by NPROCESS and round up.
    local lines_per_run=$(wc -l $GRIDLIST | \
	awk '{ x = $1/'$NPROCESS'; d = (x == int(x)) ? x : int(x)+1; print d}')

    # Use the split command to split the files into temporary files
    split --suffix-length=4 --lines $lines_per_run $GRIDLIST tmpSPLITGRID_

    # Move the temporary files into the runX-directories
    local files=$(ls tmpSPLITGRID_*)
    local i=1
    for file in $files
    do
      mv $file run$i/$GRIDLIST
      i=$((i+1))
    done
}

# Create header of progress.sh script

echo "##############################################################" > progress.sh
echo "# PROGRESS.SH" >> progress.sh
echo "# Upload current guess.log files from local nodes and check" >> progress.sh
echo "# Usage: sh progress.sh" >> progress.sh
echo >> progress.sh

# Create a run subdirectory for each process and clean up

for ((a=1; a <= NPROCESS ; a++))
do
  mkdir -p run$a
  cp $INSFILE run$a
  cd run$a ; rm -f guess.log ; rm -f $GRIDLIST ; cd ..
  echo "echo '********** Last few lines of ./run${a}/guess.log: **********'" >> progress.sh
  echo "tail ./run${a}/guess.log" >> progress.sh
done

split_gridlist

# Create PBS script to request place in queue
cat <<EOF > guess.cmd
#!/bin/bash
#PBS -l nodes=$NPROCESS
#PBS -l walltime=$WALLTIME

cd \$PBS_O_WORKDIR

cat<<EOL > startguess.sh
#!/bin/sh
cd \$PBS_O_LOCAL
mkdir \\\$((PBS_VNODENUM+1))
cd \\\$((PBS_VNODENUM+1))
cp -pr \$PBS_O_WORKDIR/run\\\$((PBS_VNODENUM+1))/* .
BINARY $INSFILE
cp -pur * \$PBS_O_WORKDIR/run\\\$((PBS_VNODENUM+1))
EOL
chmod +x startguess.sh
pbsdsh -o \$PBS_O_WORKDIR/startguess.sh
rm startguess.sh

function append_files {
    local number_of_jobs=\$1
    local file=\$2

    cp run1/\$file \$file

    local i=""
    for ((i=2; i <= number_of_jobs; i++))
    do
      cat run\$i/\$file | awk 'NR!=1 || NF==0 || \$1 == \$1+0 { print \$0 }' >> \$file
    done
}

for file in $OUTFILES
do
  append_files $NPROCESS \$file
done
EOF

# Submit job
qsub -N ${name:-"guess"} guess.cmd
