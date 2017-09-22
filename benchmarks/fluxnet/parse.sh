#! /bin/env bash

set -e
declare -A fields
fields=([TIMESTAMP]=0 [TA_ERA]=0 [P_ERA]=0 [SW_IN_ERA]=0 [NEE_CUT_MEAN]=0 [GPP_NT_VUT_REF]=0)

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
    grep -vP '^.{4}0229' "$fname" | cut -d, -f$(echo ${fields[@]} | tr ' ' ,) --output-delimiter=$'\t' |
        sed -r '1d;s/.{4}/&\t/'  > $location.csv
    rm $fname
done
