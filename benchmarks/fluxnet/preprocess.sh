#! /bin/env bash

set -e
declare -A fields
fields=([TIMESTAMP]=0 [TA_ERA]=0 [SW_IN_ERA]=0 [P_ERA]=0 [NEE_CUT_MEAN]=0 [GPP_NT_VUT_REF]=0)

iconv -f latin1 -t utf-8 < *_sitelist.csv |
    awk -vFPAT='[^,]*|"[^"]*"' -vOFS='\t' '$1 !~ /^$/ { print $7,$6,$1}' > fixed.txt

for f in *.zip; do
    fname=$(basename "$f")
    location=${fname:4:6}
    unzip -q $f '*_FULLSET_DD_*'
    fname=$(ls *_DD_*.csv)
    i=1
    for field in $(head -1 *_DD_*.csv | tr , ' '); do
        if [ ${fields[$field]+_} ]; then
            fields[$field]=$i
        fi
        i=$(( $i+1 ))
    done
    grep $location fixed.txt | tee $location.csv |
        awk -vOFS='\t' 'function norm(x, y){y=int(x*2); return(2*x<y?y-1:y)/2+.25}
            {print norm($1), norm($2), $3}' >> gridlist.txt
    grep -vP '^.{4}0229' "$fname" | cut -d, -f$(echo ${fields[@]} | tr ' ' ,) --output-delimiter=$'\t' |
        sed -r '1d;s/.{4}/&\t/' > tmp
    cut -f-5 tmp >> $location.csv
    cut -f6- tmp | sed "s/^/$location\t/" >> validation.csv
    rm $fname
done
rm tmp fixed.txt
