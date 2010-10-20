# This file contains a function for creating a scatter plot using gnuplot.
#
# Parameters:
# $1 title of the graph
# $2 title of the x axis
# $3 title of the y axis
# $4 the input data file with two columns
# $5 the output filename
function scatter_plot {
    local title=$1
    local xtitle=$2
    local ytitle=$3
    local datafile=$4
    local image=$5

# The gnuplot script
    cat <<EOF | gnuplot &> /dev/null
set title "${title}"
set xlabel "${xtitle}"
set ylabel "${ytitle}"
set key off
set grid ytics
set terminal png
set output "${image}"
set fit logfile "/dev/null"
y(x) = m*x+c
fit y(x) "${datafile}" via m,c
set label 'y = %.3f',m,'x + %.3f',c at screen 0.15,0.85 front
plot "${datafile}" with points lt 3 pt 5, y(x)
EOF
}
