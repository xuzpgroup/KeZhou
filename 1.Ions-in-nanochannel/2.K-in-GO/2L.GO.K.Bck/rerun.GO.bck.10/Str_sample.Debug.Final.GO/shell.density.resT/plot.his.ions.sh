#!/bin/bash
## plot figure 
cat << EOF > plot.plt
#set terminal postscript eps color solid linewidth 2 "Helvetica" 26
# set terminal pngcairo  transparent enhanced font "arial,12" fontscale 1.0  ##size 512, 280 
set term png  enhanced  font "arial,15" fontscale 1.0  ##size 512, 280 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'His.shell.z.png'
set style fill solid
set multiplot layout 2,1 title ""


#set xrange [0:1000]
#set yrange [-2:6]
#set yrange [0:5]
set xlabel 'Position Z({\305})'
#set ylabel 'Q ({/Symbol-Oblique w})'
set ylabel 'Distribution probability'
set xtics  2
set ytics  0.2
#unset xlabel
#unset xtics
#set key center at 130, 4.0  ## set legend position
#set logscale

#p 'Q.omega.data' u (\$1*2*${PI}/${T}*1.0E-3):(\$3) w p pt 0 lw 2 lc 0 lt 1  title ''
p './1_Distribution.GO.10.0/Water.distribution.data' u (\$1):(\$2) w l lw 2 lc 1 lt 1  title 'O atom, water',\
'./His.shell.z.dat' u (\$1):(\$2) w l lw 2 lc 2 lt 1  title 'O atom, shell',\
'./His.shell.z.dat' u (\$1):(\$4) w l lw 2 lc 3 lt 1  title 'Mass',\





EOF
gnuplot plot.plt
rm plot.plt


