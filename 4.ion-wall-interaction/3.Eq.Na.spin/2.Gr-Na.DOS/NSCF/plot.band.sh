#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set terminal postscript eps enhanced color  font 'Arial,20' linewidth 1 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'band.NaG.eps'
set style fill solid
set multiplot layout 1,2 title ""
set size ratio 1.8

load "mycol.pal"

# set lmargin 8 
# set rmargin 2 
# set key font ',20'
# set key b

# set xlabel '{/Arial-Italic E} (eV)'
# set ytics  2
# set xtics  4
# set xrange [2.0:8.0]
set yrange [:15]
set ytics  10
# set xtics  1
set ylabel '{/Arial-Italic E} (eV)'
unset xlabel
unset xtics
set xtics ("G" 0.0, "M"  0.36738021, "K" 0.57948725, "G" 1.0037013)


set key l
# p 'DOS0' u (\$1*1):(\$2)  w l ls 21 lw 4  title '',\

# p 'BAND.dat' u (\$1*1):(\$2)  w lp ls 21 pt 7 lw 4  title '',\

p 'BAND.dat' u (\$1*1):(\$2)  w l  ls 21  lw 4  title 'up',\


p 'BAND.dat' u (\$1*1):(\$3)  w l  ls 21  lw 4  title 'down',\








EOF
gnuplot plot.plt
rm plot.plt

file="DOS.graphene"
# convert -density 500 ${file}.eps ${file}.png 








