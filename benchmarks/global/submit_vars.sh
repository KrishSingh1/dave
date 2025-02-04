NPROCESS=24
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=96
elif (( "$ARCH" == "levante" )); then
    NPROCESS=256
fi
