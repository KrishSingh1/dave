NPROCESS=1
if (( "$ARCH" == "levante" )); then
    PARTITION=shared
    NPROCESS=1
    NTASKSPERNODE=1
    WALLTIME="1:00:00"
    WALLTIME_APPEND="0:30:00"
fi
