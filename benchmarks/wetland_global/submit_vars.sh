NPROCESS=24
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=96
elif [[ "$ARCH" == "levante" ]]; then
    NPROCESS=128
    WALLTIME=02:00:00 # took 01:13:21
    WALLTIME_APPEND=0:30:00 # took 0:07:31
 fi
