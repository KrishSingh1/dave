#! /bin/env bash

set -e
fields=(TIMESTAMP TA_ERA SW_IN_ERA P_ERA NEE_VUT_REF GPP_NT_VUT_REF)

function indexes {
    local fname=$1
    shift
    declare -A names
    for n in "$@"; do
        names[$n]=0
    done

    local i=1
    for field in $(head -1 $fname | tr , ' '); do
        if [ ${names[$field]+_} ]; then
            names[$field]=$i
        fi
        i=$(( $i+1 ))
    done
    echo "${names[@]}" | tr ' ' ,
}

iconv -f latin1 -t utf-8 < *_sitelist.csv |
    awk -vFPAT='[^,]*|"[^"]*"' -vOFS='\t' '$1 !~ /^$/ {print $7,$6,$1}' > fixed.txt

for f in raw/*.zip; do
    fname=$(basename "$f")
    location=${fname:4:6}
    if [ -f $location.csv ]; then
        echo skipping $fname - duplicate at $location
        continue
    fi
    unzip -q $f '*_FULLSET_DD_*'  '*_FULLSET_MM_*'
    if ! grep -q ,NEE_VUT_REF, <(head -1 *_MM_*.csv); then
        echo skipping $location - missing NEE
        rm *_{DD,MM}_*.csv
        continue
    fi
    grep -vP '^.{4}0229' *_DD_*.csv |
        cut -d, -f$(indexes *_DD_*.csv ${fields[@]}) --output-delimiter=$'\t' |
        sed -r '1d;s/.{4}/&\t/' > tmp

    if (( $(wc -l tmp | cut -f1 -d' ') % 365 != 0 )); then
        echo "skipping $location"
    else
        grep $location fixed.txt | tee $location.csv |
            awk -vOFS='\t' 'function norm(x, y){y=int(x*2); return(2*x<y?y-1:y)/2+.25}
                {print norm($1), norm($2), $3}' >> gridlist.txt
        cut -f-5 tmp >> $location.csv
        cut -f6- tmp | sed "s/^/$location\t/" >> daily.csv

        cut -d, -f$(indexes *_MM_*.csv ${fields[@]: -3:3}) --output-delimiter=$'\t' *_MM_*.csv |
            sed -r "1d;s/^/$location\t/" >> monthly.csv
    fi
    rm *_{DD,MM}_*.csv
done
rm tmp fixed.txt
