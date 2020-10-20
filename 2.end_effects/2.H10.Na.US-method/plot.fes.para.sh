#!/bin/bash
## plot figure 
KJ2eV=0.0103643
MidP=24.75
cat << EOF > plot.plt
# set term post enhanced 20 color eps
# set encoding iso_8859_1
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'fes.para.png'
set style fill solid
set multiplot layout 1,1 title ""


#set xlabel '{/Times-Italic d} [{\305}]'
set xlabel 'd [{\305}]'
# set ytics  2
# set xtics  4
# set xrange [6.0:24.0]
# set yrange [-0.9:6]
#set ytics  5
# set xtics  1
#set ylabel '{/Times-Italic F} [kJ/mol]'
#set ylabel '{/Times-Italic F} [eV]'
set ylabel 'F [KJ/mol]'
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
#set logscale y 
# unset label
# set label 1 at 7, 5.2 'Bulk SPC/E, our = 2.75'
# set label 2 at 7, 5.5 'Bulk TIP3P, PCCP = 6.22'
# set arrow  from 6.0, 2.75 to 24.0, 2.75 nohead lt 2 lc rgb '#DC143C'  ## ad a line
set key b
p './His/free.out' u (\$1*10):(\$2) w l lw 4 lc 1 lt 1  title '5ns, pre 250 ps',\
'./His.frq.2/free.out' u (\$1*10):(\$2) w l lw 4 lc 2 lt 1  title '5ns, pre 500 ps',\
'./His.frq.4/free.out' u (\$1*10):(\$2) w l lw 4 lc 3 lt 1  title '5ns, pre 1000 ps',\
'./His.3ns/free.out' u (\$1*10):(\$2) w l lw 4 lc 4 lt 1  title '3ns, pre 250 ps',\
'./His.windos/free.out' u (\$1*10):(\$2) w l lw 4 lc 5 lt 1  title 'half windows',\



### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt










