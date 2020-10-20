#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced  font "arial,15" fontscale 1.0  ##size 512, 280 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'msd.all.cases.GO.1L.K.png'
set style fill solid
set multiplot layout 1,1 title "K^+ in 1L/GO"

#set format x "%.e "
#set yrange [0:3.2]
set xlabel 't [ps]'
#set ytics  5
set xtics  500
set xrange [0:1000]
set yrange [0:5]
set ytics  2.5
set ylabel 'MSD [{\305}^2]'
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
p '2_MSD_for_gas/msd.data' u (\$1):(\$2) w lp pt 0 lw 2 lc 1 lt 1  title 'case 1',\
'rerun.5/2_MSD_for_gas/msd.data' u (\$1):(\$2) w lp pt 0 lw 2 lc 2 lt 1  title 'case 2',\
'rerun.10/2_MSD_for_gas/msd.data' u (\$1):(\$2) w lp pt 0 lw 2 lc 3 lt 1  title 'case 3',\
'rerun.15/2_MSD_for_gas/msd.data' u (\$1):(\$2) w lp pt 0 lw 2 lc 4 lt 1  title 'case 4',\
'rerun.20/2_MSD_for_gas/msd.data' u (\$1):(\$2) w lp pt 0 lw 2 lc 5 lt 1  title 'case 5',\
'rerun.25/2_MSD_for_gas/msd.data' u (\$1):(\$2) w lp pt 0 lw 2 lc 6 lt 1  title 'case 6',\
'rerun.30/2_MSD_for_gas/msd.data' u (\$1):(\$2) w lp pt 0 lw 2 lc 7 lt 1  title 'case 7',\
'rerun.35/2_MSD_for_gas/msd.data' u (\$1):(\$2) w lp pt 0 lw 2 lc 8 lt 1  title 'case 8',\
'rerun.40/2_MSD_for_gas/msd.data' u (\$1):(\$2) w lp pt 0 lw 2 lc 9 lt 1  title 'case 9',\

#'rerun.45/2_MSD_for_gas/msd.data' u (\$1):(\$2) w lp pt 0 lw 2 lc 10 lt 1  title 'case 10',\
#'msd.ave.dat' u (\$1):(\$2) w lp pt 0 lw 4 lc 0 lt 1  title 'ave',\

 

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt


