#!/bin/bash
## plot figure 
dH=1.0
cat << EOF > plot.plt
set terminal postscript eps enhanced color  font 'Arial,20' linewidth 1 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'fat.band.eps'
set style fill solid
set multiplot layout 1,2 title ""
set size ratio 1.8

load "mycol.pal"

set lmargin 8 
set rmargin 2 
set key font ',20'  box opaque
# set key at 0.15,-17.0
# set key b

# set xlabel '{/Arial-Italic E} (eV)'
# set ytics  2
# set xtics  4
# set xrange [2.0:8.0]
set yrange [-30:15]
# set ytics  2
# set xtics  1
set ylabel '{/Arial-Italic E} (eV)'
unset xlabel
unset xtics
# set xtics ("G" 0.0, "M"  1.46952088, "K" 2.31794907, "G" 4.01480545)
set xtics ("G" 0.0, "M"  0.36738021, "K" 0.57948725, "G" 1.0037013)

set style fill transparent solid 0.8 noborder

set key l
set key b
# p 'DOS0' u (\$1*1):(\$2)  w l ls 21 lw 4  title '',\

# p 'BAND.dat' u (\$1*1):(\$2)  w lp ls 21 pt 7 lw 4  title '',\




p 'PBAND_Na_DW.dat' u (\$1*1):(\$2)  w l  ls 21  lw 1  title 'up',\
'PBAND_Na_DW.dat'  u (\$1):(\$2):(\$3*5) \
w points ls 22 pt 7 ps var title 'Na ,s' ,\

p 'PBAND_Na_UP.dat' u (\$1*1):(\$2)  w l  ls 21  lw 1  title 'down',\
'PBAND_Na_UP.dat'  u (\$1):(\$2):(\$3*5) \
w points ls 22 pt 7 ps var title 'Na, s' ,\









EOF
gnuplot plot.plt
rm plot.plt

file="fat.band"
convert -density 500 ${file}.eps ${file}.pts.SC.NaG.png 








