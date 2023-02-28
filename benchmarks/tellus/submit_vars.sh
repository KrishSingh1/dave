# the tellus bm run on GAP keal's 8 sandy nodes each with 16 CPUs, took approx 25 hours
NPROCESS=30
if [[ $ARCH == "aurora" ]]; then
    NPROCESS=100
fi
