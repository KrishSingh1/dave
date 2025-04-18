NPROCESS=6
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=20
elif [[ "$ARCH" == "levante" ]]; then
    PARTITION=shared
    NPROCESS=20
    NTASKSPERNODE=20
    WALLTIME="1:00:00"
    WALLTIME_APPEND="0:30:00"
fi
