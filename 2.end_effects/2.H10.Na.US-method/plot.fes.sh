#!/bin/bash
## plot figure 
KJ2eV=0.0103643
MidP=24.75
KJmol2kbT=2.493
cat << EOF > plot.plt
# set term post enhanced 20 color eps
# set encoding iso_8859_1
# set term png  enhanced
set terminal postscript eps enhanced color  font 'Arial,25' linewidth 1 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'FES.H10.M1.Na.eps'
set style fill solid
set multiplot layout 1,1 title ""
set size ratio 0.72


#set xlabel '{/Times-Italic d} [{\305}]'
set xlabel '{/Arial-Italic x} (nm)'
set ytics  1
set xtics  0.5
set xrange [-1.0:1.0]
set yrange [-0.3:2]
#set ytics  5
# set xtics  1
#set ylabel '{/Times-Italic F} [kJ/mol]'
#set ylabel '{/Times-Italic F} [eV]'
# set ylabel 'F [KJ/mol]'
set ylabel 'PMF ( {/Arial-Italic k}_B{/Arial-Italic T} )'
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
#set logscale y 
# unset label
# set label 1 at 7, 5.2 'Bulk SPC/E, our = 2.75'
# set label 2 at 7, 5.5 'Bulk TIP3P, PCCP = 6.22'
# set arrow  from 6.0, 2.75 to 24.0, 2.75 nohead lt 2 lc rgb '#DC143C'  ## ad a line

set style fill transparent solid 0.4 noborder
set y2range [0:8.0]
set xrange [-1.0:1.0]
set x2range [-1.0:1.0]
set y2label '{/Arial-Italic N}_c'
set y2tics 2
set ytics nomirror

p 'free.out' u (\$1):(\$2/${KJmol2kbT}) \
w lp pt 7 ps 1 lw 4 lc rgb 'red'  title 'Na^+, PMF',\
'./Nc.path.dat' u (\$1/10-3.15):(\$2) axes x2y2  \
w p pt 5 lc rgb 'black'  title '{/Arial-Italic N}_c',\

# '../1.H10.Na.2/Distribution.width.data/Nc.path.dat' \
# u (\$1/10.0-3.15):(\$2) axes x2y2 \
# w p pt 5 lc rgb 'black'  title '{/Arial-Italic N}_c',\

# 'free.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 7 lw 8 lc 2 lt 1  title '' smooth  sbezier


# p 'log' u (\$1*10):(\$2) w l lw 1 lc 1 lt 1  title ''

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt










