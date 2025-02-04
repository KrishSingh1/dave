NPROCESS=6
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=1
elif (( "$ARCH" == "levante" )); then
    PARTITION="shared"
    WALLTIME=01:00:00
    NPROCESS=7
fi
