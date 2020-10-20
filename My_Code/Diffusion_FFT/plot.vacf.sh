#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

set output 'VACF.png'
set style fill solid
set multiplot layout 1,1 title ""

set xrange [0.0:10]
set xtics 1
set ytics 0.5
set grid

#set yrange [1:3.9]
#set yrange [0:200]
set xlabel 'Time (ps)'
set ylabel 'VACF'
p '1nor_facf.data' u (\$1*1.0E-2):(\$2) w lp pt 0 lw 2 lc 1 lt 1  title ''

EOF
gnuplot plot.plt
rm plot.plt





