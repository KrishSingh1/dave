#!/bin/bash

dir=$(grep flux_dir common/../paths.ins | awk -vFPAT='"[^"]*"' 'END{gsub(/"/, "", $2); print $2}')

source scatter_plot.sh

function compare {
    local var=$1
    local outfile=fluxnet_${var,,}.png
    local mod_idx=$(( $2 + 3 ))

    paste <(cut -f$2 ${dir}validation.csv) <(awk "NR > 1 {print \$$mod_idx}" dfluxnet.out) |
        grep -v "\-9999" > tmp

    scatter_plot $var Observed Modelled tmp $outfile

    rm tmp
    describe_image $outfile "Daily $var: FLUXNET vs. LPJ-GUESS. Units: kgC m-2." embed
}

describe_benchmark "LPJ-GUESS - FLUXNET Benchmarks (global PFTs)"

compare NEE 2
compare GPP 3
