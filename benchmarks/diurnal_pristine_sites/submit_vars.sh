NPROCESS=5
INPUT_MODULE=watch_diurnal
if (( "$ARCH" == "levante" )); then
    PARTITION=shared
    NPROCESS=5
    NTASKSPERNODE=5
    WALLTIME="1:00:00"
    WALLTIME_APPEND="0:10:00"
fi
