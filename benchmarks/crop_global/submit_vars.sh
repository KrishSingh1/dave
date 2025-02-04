NPROCESS=30
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=96
elif [[ "$ARCH" == "levante" ]]; then
    NPROCESS=256
    WALLTIME=8:00:00
fi
