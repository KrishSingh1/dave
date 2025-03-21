#!/bin/bash

# Subroutine for deltareport.sh
# By Johan Nord, 2017.
# Comparison of data in the benchmark report catalog with comparator benchmark.
# Run in the BM outputfolder where there is a report directory to be compared.
# Output to (new) catalog report_delta.

# Input parameters:
#     $1 = $diffreportpath, i.e. path to where to put the comparison report are being generated. It is normally sitting in the pwd of the call to the deltareport script.
#          (The path to where the report to compare with resides, is grabbed from the symlink in the catalog on the $diffreportpath)

# N.B.:
#     If there is no difference beween report gmap and reference-comparator gmap, then no delta-gmap is made.


set -e

# Set paths to work on
reportpath=report
diffreportpath=$1			# e.g. report_delta or a permutation of it. in v1&2 it was the folder name. Now in v3 it is the full path.
comparatorpath=`cd ${diffreportpath}/comparator_benchmark; cd $(pwd -P); pwd -P`       #do not change string constant 'comparator_benchmark' unless changed in deltareportv1.sh
benchmark=$(basename $(pwd -P))
postprocesssh="common/../${benchmark}/postprocess.sh"

# Make paths absolute
reportpath=`cd $reportpath; pwd -P`
diffreportpath=`cd $diffreportpath; pwd -P`

# Diagnostics
echo PWD: $(pwd -P)
echo -e reportpath '\t' `ls $reportpath -d`
echo -e comparatorpath '\t' `ls $comparatorpath -d`
echo -e diffreportpath '\t' `ls $diffreportpath -d`
echo -e postprocesssh '\t' `ls $postprocesssh`

set +e


### This section deals with making delta images for the gmapall- and gmap section of benchmarks xml file

# Make a list of the relevant gmapall-tslices based on the roots above (or do this list as part of the job above)
   gmapall_files=$(grep gmapall  $postprocesssh | cut -d\  -f2)                   # i.e. the tslices of gmapall calls in the BM's postprocess.sh  e.g. cmass1961to1990.txt
gmapall_prefixes=$(grep gmapall  $postprocesssh | cut -d\  -f4)                # e.g. cmass_
 gmapall_legends=$(grep gmapall  $postprocesssh | cut -d\  -f6)
 gmapall_n_files=$(grep gmapall  $postprocesssh | wc -l)
echo do delta and gmapall on $gmapall_files				                        # e.g. cmass1961to1990.txt
       echo gmapall_prefixes $gmapall_prefixes
        echo gmapall_legends $gmapall_legends

# Make a list of the relevant gmap-tslices (NOT gmapall) based on the roots above
   gmapNONall_files=$(grep "gmap "  $postprocesssh | cut -d\  -f2) 			# Removed | grep "1990.txt" --- Correct? Påverkan på andra benchmarks=?
gmapNONall_outfiles=$(grep "gmap "  $postprocesssh | sed 's|-o |\n-o |g' | grep ^\-o | cut -d\  -f2)
  gmapNONall_titles=$(grep "gmap "  $postprocesssh | sed 's|-t |\n-t |g' | grep ^\-t | sed 's|\ |_|g' | cut -d\' -f2)
 gmapNONall_n_files=$(grep "gmap "  $postprocesssh | wc -l)			# Removed | grep "1990.txt" --- Correct?

echo and do delta and gmap on $gmapNONall_files			# e.g. cmass1961to1990.txt
echo files  $gmapNONall_files
echo titles $gmapNONall_titles
echo outfiles  $gmapNONall_outfiles
echo n files  $gmapNONall_n_files


# Do delta on tslice 1961to1990,
# incl check column headers are identical (todo?)
# ideally these resulting delta tslices should reside in the $diffreportpath

echo -n to file.delta:\
files_to_delta="$gmapall_files $gmapNONall_files"
for file in $files_to_delta ; do

  if [ ! "$file" == " " ]; then
    file1=${reportpath}/../$file
    file2=${comparatorpath}/../$file
    outfile=${diffreportpath}/$file.delta

#    if diff -q file1 file2 1>/dev/null; then
#      echo "identical" > $outfile
#    else
      echo -n ${file}\
      #echo delta 1 ${file1}   2 ${file2}  - to - ${outfile}

      delta ${file1} ${file2} -i Lon -i Lat -o ${outfile}  1>/dev/null
#    fi
  fi
done
echo

# gmapalls of tslices of the relevant-tslices list. This generates delta imagess in $diffreportpath
echo gmapalls of tslices of the relevant-tslices list.
pwd
echo $diffreportpath

cd $diffreportpath

for i in `seq $gmapall_n_files`; do

  file=$(echo $gmapall_files | cut -d\  -f$i)
  prefix=$(echo $gmapall_prefixes | cut -d\  -f$i)
  legend=$(echo $gmapall_legends | cut -d\  -f$i)

  echo -n ${file}.delta..$prefix\
  gmapall $file.delta -P $prefix -portrait -slog -c GREYRED DEPTH	1>/dev/null
  # detta skall funka! även om vill i framtiden ändra DEPTH till en färgskala som har DEPTHs överstvärde=grön
  # och kanske tydligare skillnad i färgerna allra närmast noll (varav ev gul på negativa sidan))

done
echo

cd ..

# gmaps (non-gmapalls) of tslices of the relevant-tslices list. This generates delta imagess in $diffreportpath
echo  gmaps of tslices of the relevant-tslices list.

cd $diffreportpath
echo filesnames: $gmapNONall_files
echo titles: $gmapNONall_titles
echo outfilename: $gmapNONall_outfiles

for i in `seq $gmapNONall_n_files`; do

  #echo i = .$i.
  #i = $(tr -d \  $i)
  #echo i = .$i.
  file=$(echo $gmapNONall_files | cut -d\  -f$i)
  title=$(echo $gmapNONall_titles | cut -d\  -f$i)
  #  legend=$(echo $gmapNONall_legends | cut -d\  -f$i)
  outfile=$(echo $gmapNONall_outfiles | cut -d\  -f$i)


  # Make gmaps of deltas, or in a few special cases, make diff-detection images
  #echo -n $i ${file}.delta

  #if [ $file == "biomes_lai1961to1990.txt" ] || [ $file == "lai1961to1990max.txt" ]; then
    #compare ${reportpath}/biomes.jpg ${comparatorpath}/biomes.jpg $outfile   1>/dev/null
    #compare ${reportpath}/$outfile ${comparatorpath}/$outfile $outfile   1>/dev/null
  #else
    gmap $file.delta -t $title -lon 1 -lat 2 -i 3 -portrait -o $outfile -pixoffset 0.0 0.0 -slog -c GREYRED DEPTH   1>/dev/null
  #fi

  if [ ! -f $outfilec ]; then
    convert -background white -size 1500x -gravity Center -weight 700 -pointsize 200 caption:"No difference. (Or script bug)" $outfile
  fi

done
echo

cd ..


### This section deals with making thumbnails for all images

# Make thunbmnails of all images in $reportpath, $comparatorpath, $diffreportpath and put them in $diffreportpath/{tmbs,tmbs_cmp,tmbs_delta} 
# AND edit the $diffreportpath/report.xml markup as to replace textlinks to the images with thubnails as links
# that point to images in $reportpath, $comparatorpath, $diffreportpath.
# This is the trickiest part because it requires to know how to markup this in the report.xml format.
# N.B. It is Not needed to edit the file $diffreportpath/report.xml: the correct img src url is already in the report.xml file


# Make thumbnails of ALL jpgs in the report_delta folder

allimages=$(grep 'image src=' ${reportpath}/report.xml | cut -d\" -f2)
plots=$(grep 'image src=' ${reportpath}/report.xml | grep 'embed="1"' | cut -d\" -f2)	# Line contains 'embed=', assume that it is a plot

echo all images: $allimages
echo plots: $plots

cd $diffreportpath
mkdir -p thmbs_delta
mkdir -p thmbs_org
mkdir -p thmbs_ref

for imgfile in $allimages; do

 if [ "$benchmark" = "europe" ] || [ "$benchmark" = "pristine_sites" ] || [ "$benchmark" = "diurnal_pristine_sites" ] || [ "$benchmark" =  "emdi_europe" ] ; then
    convert -verbose ${reportpath}/$imgfile -resize 15% thmbs_org/${imgfile}.thmb.jpg      1>/dev/null
    convert -verbose ${comparatorpath}/$imgfile -resize 15% thmbs_ref/${imgfile}.thmb.jpg  1>/dev/null
    convert -verbose $imgfile -resize 15% thmbs_delta/${imgfile}.thmb.jpg   1>/dev/null

 else
    convert -verbose ${reportpath}/$imgfile -resize 15% -crop 100x50%+0+120% thmbs_org/${imgfile}.thmb.jpg   1>/dev/null
    convert -verbose ${comparatorpath}/$imgfile -resize 15% -crop 100x50%+0+120% thmbs_ref/${imgfile}.thmb.jpg  1>/dev/null
    convert -verbose $imgfile -resize 15% -crop 100x50%+0+120% thmbs_delta/${imgfile}.thmb.jpg   1>/dev/null

 fi

done

echo $allimages

cd ..

### This section deals with adjusting the index.htlm file
# Edit the index.html

cd $diffreportpath

 # gör alla >div> till ny rad och formatera fint: grep tar bort blankrader. Sista sed för att separera svn.info från tabellen ovan.
cat index.html | sed 's|<div|\n<div|g' - | sed 's|</div>|</div\>\n|g' | sed 's|<h3>Meta data</h3>|\n<h3>Meta data</h3>|g' - | grep -v -e '^$' >index.html.formatted.tmp

 # sen skapa denna raden med ursprungsraden som mall:

echo \  >index.html
i=1			# Counter for the jpeg files to be looked for. We look only once for each jpeg file.
while IFS= read -r line; do

  imgfile=$(echo $allimages | cut -d\  -f$i)

  #echo $i "$imgfile"		# For debugging

  if [[ $line == *"$imgfile"* ]] && [ -n "$imgfile" ]; then		# Line contains image file => replace line

   if [[ $plots == *"$imgfile"* ]] && [ -n "$imgfile" ]; then		# Line contains plot file => do not show thumbnail

    #echo $line			# For debugging

    lastpartoftheline=$(echo $line | cut -d\> -f4-)
    echo "<div><image src=\"../report/$imgfile\"/><image src=\"comparator_benchmark/$imgfile\"/><br/>${lastpartoftheline}" >> index.html

   else									# Line probaly contains gmap jpeg => show thumbnail

    echo "<div>$imgfile<br> <a href=\"../report/$imgfile\"><IMG HEIGHT=270 WIDTH=400 SRC=\"thmbs_org/$imgfile.thmb.jpg\"</a><a href=\"comparator_benchmark/$imgfile\"><IMG HEIGHT=270 WIDTH=400 SRC=\"thmbs_ref/$imgfile.thmb.jpg\"</a><a href=\"$imgfile\"><IMG HEIGHT=270 WIDTH=400 SRC=\"thmbs_delta/$imgfile.thmb.jpg\"</a></div>" >>index.html   

   fi

   let "i++"

  else									# Line does not contain jpeg => copy line from file

    echo "$line" >> index.html

  fi

done < "index.html.formatted.tmp"

rm index.html.formatted.tmp


sed -i "s|Comparator benchmark|\<p\>Comparator benchmark|" index.html

echo Finished OK
