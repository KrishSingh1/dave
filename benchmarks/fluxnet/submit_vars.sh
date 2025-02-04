NPROCESS=5
INPUT_MODULE=fluxnet
if (( "$ARCH" == "levante" )); then
    PARTITION=shared
    NPROCESS=20
    NTASKSPERNODE=20
    WALLTIME="0:10:00"
    WALLTIME_APPEND="0:10:00"
fi
