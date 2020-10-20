#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

set output 'Trap.png'
set style fill solid
set multiplot layout 1,1 title ""

#set xrange [0.0:4]
#set xtics 1
#set ytics 0.5
set grid
#set logscale y
set yrange [-1:1]
set xrange [0:600]
set xlabel 'Time (ps)'
set ylabel 'D(10^-^1^1m^2/s)'

p '3tc.data' u (\$1*1.0E3):(\$2) w lp pt 0 lw 2 lc 1 lt 1  title ''

EOF
gnuplot plot.plt
rm plot.plt








