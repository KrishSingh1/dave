NPROCESS=12
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=48
elif [[ "$ARCH" == "levante" ]]; then
    NPROCESS=128
    WALLTIME=03:00:00
fi

