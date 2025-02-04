NPROCESS=24
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=96
elif (( "$ARCH" == "levante" )); then
    NPROCESS=128
    WALLTIME=02:00:00
 fi
