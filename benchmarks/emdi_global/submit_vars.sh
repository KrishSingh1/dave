NPROCESS=6
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=48
elif [[ "$ARCH" == "levante" ]]; then
    PARTITION=shared
    NPROCESS=20
    NTASKSPERNODE=20
    WALLTIME="00:30:00"
    WALLTIME_APPEND="0:02:00"
fi

