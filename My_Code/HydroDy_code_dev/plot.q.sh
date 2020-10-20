#!/bin/bash
## plot figure 
KJ2eV=0.0103643
KJmol2kbT=2.493
MidP=24.75
cat << EOF > plot.plt
# set term post enhanced 20 color eps
# set encoding iso_8859_1
# set term png  enhanced 
set terminal postscript eps enhanced color  font 'Arial,32' linewidth 1
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Charge.profile.AIMD.eps'
set style fill solid
set multiplot layout 1,1 title ""
set size ratio 0.4

#set xlabel '{/Times-Italic d} [{\305}]'
set xlabel '{/Arial-Italic z} (nm)'
set xtics  0.20
# set ytics  2
set mxtics 2
set mytics 2
set xrange [-0.5:1.5]
set xrange [0.1:0.9]
# set yrange [:6]

# unset ylabel
# set lmargin 2
# set rmargin 2
# set ylabel 'F [KJ/mol]'
set ylabel 'Q (e)'
#set ytics  0.4
#unset xlabel
# unset ytics
#set key center at 0.55, 0.8  ## set legend position
#set logscale y 
# unset label
set key font ",25"
# set key l

set style fill solid 0.5

# set label 1 at 7, 5.2 'Bulk SPC/E, our = 2.75'
# set label 2 at 7, 5.5 'Bulk TIP3P, PCCP = 6.22'
# set arrow  from 6.0, 2.75 to 24.0, 2.75 nohead lt 2 lc rgb '#DC143C'  ## ad a line

#p 'fes.dat' u (\$1*10):(\$2/${KJmol2kbT}*0.0103643) w lp pt 4 ps 2 lw 2 lc 1 lt 1  title 'd',\

p 'charge.density.profile.dat' u (\$1*.1-0.4):(\$2) \
 w filledcurves lw 4 lc rgb '#006699' lt 1  title '',\

 




### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt


# for i in *.eps 
# do \
echo "Convert to png"
# convert -density 500 Density.profile.AIMD.eps Density.profile.AIMD.png 
# done 





