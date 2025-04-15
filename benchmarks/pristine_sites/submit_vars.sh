NPROCESS=5
if (( "$ARCH" == "levante" )); then
    PARTITION=compute
    NPROCESS=5
    NTASKSPERNODE=5
    WALLTIME="0:10:00"
    WALLTIME_APPEND="0:02:00"
fi


