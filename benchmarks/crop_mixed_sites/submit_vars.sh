NPROCESS=6
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=20
elif (( "$ARCH" == "levante" )); then
    PARTITION=shared
    NPROCESS=20
    NTASKSPERNODE=20
    WALLTIME="0:10:00"
    WALLTIME_APPEND="0:02:00"
fi
