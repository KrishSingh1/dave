NPROCESS=12
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=48
elif [[ "$ARCH" == "levante" ]]; then
    NPROCESS=128
    WALLTIME=03:00:00 # took 01:08:56 on Levante 
    WALLTIME_APPEND="0:10:00" # took 0:02:26 on Levante
fi
