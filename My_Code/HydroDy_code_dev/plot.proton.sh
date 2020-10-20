#!/bin/bash
## plot figure 
cat << EOF > plot.plt
# set term post enhanced 20 color eps
# set encoding iso_8859_1
# set term png  enhanced 
set terminal postscript eps enhanced color  font 'Arial,16' linewidth 1
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Proton.AIMD.eps'
set style fill solid
set multiplot layout 1,1 title ""
set size ratio 0.2

#set xlabel '{/Times-Italic d} [{\305}]'
set xlabel '{/Arial-Italic t} (ps)'
set ytics  2
# set mxtics 2
# set mytics 2
set yrange [:6]

# unset ylabel
# set lmargin 2
# set rmargin 2
# set ylabel 'F [KJ/mol]'
#set ytics  0.4
#unset xlabel
# unset ytics
#set key center at 0.55, 0.8  ## set legend position
#set logscale y 
# unset label
# set key font ",25"
# set key l

set style fill solid 0.5

# set label 1 at 7, 5.2 'Bulk SPC/E, our = 2.75'
# set label 2 at 7, 5.5 'Bulk TIP3P, PCCP = 6.22'
# set arrow  from 6.0, 2.75 to 24.0, 2.75 nohead lt 2 lc rgb '#DC143C'  ## ad a line

#p 'fes.dat' u (\$1*10):(\$2/${KJmol2kbT}*0.0103643) w lp pt 4 ps 2 lw 2 lc 1 lt 1  title 'd',\

set ylabel '{/Arial-Italic n}_{H^+}'
# set xrange [0.9:1.2]
# set xtics  0.1
# set ytics  0.04

p 'proton_sequence.dat' u (\$1*0.5E-3):(\$2) \
 w l lw 2 lc rgb 'blue' lt 1  title '',\
 




### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt


# for i in *.eps 
# do \
echo "Convert to png"
# convert -density 500 Density.profile.AIMD.eps Density.profile.AIMD.png 
# done 





