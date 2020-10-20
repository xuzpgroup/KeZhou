#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'gas.tar.png'
set style fill solid
set multiplot layout 1,1 title ""


#set ytics 0.1
#set xrange [0.2:1.0]
#set xrange [0.0:1.0]
#set yrange [1:3.9]
#set yrange [0:200]
#set xrange [0:1000]
#set yrange [0:3.2]
set xlabel 'X [{\305}]'
set ytics  5
set xtics  5
set ylabel 'Y [{\305}]'
set size square
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
p 'gas.tar.dat' u (\$2):(\$3) w lp pt 4 lw 2 lc 3 lt 1  title '' ,\



### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt

