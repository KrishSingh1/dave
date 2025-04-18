NPROCESS=1
if [[ "$ARCH" == "levante" ]]; then
    PARTITION=shared
    NPROCESS=1
    NTASKSPERNODE=1
    WALLTIME="0:10:00"
    WALLTIME_APPEND="0:02:00"
fi
