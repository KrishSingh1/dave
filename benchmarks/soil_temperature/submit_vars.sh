NPROCESS=6
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=48
elif [[ "$ARCH" == "levante" ]]; then
    NPROCESS=128
    WALLTIME=00:10:00 # takes only 2:30 on Levante
fi
