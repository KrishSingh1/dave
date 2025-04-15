NPROCESS=6
if [[ "$ARCH" == "lunarc" ]]; then
    NPROCESS=1
elif (( "$ARCH" == "levante" )); then
    PARTITION="shared"
    WALLTIME=00:10:00 # takes only 39 seconds on Levante
    WALLTIME_APPEND="0:02:00"
    NPROCESS=7    
fi
