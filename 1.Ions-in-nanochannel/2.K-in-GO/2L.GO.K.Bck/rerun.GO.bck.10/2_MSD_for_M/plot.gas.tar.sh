X=63.000000
Y=111.544072

#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced  font "arial,15" fontscale 1.0  ##size 512, 280 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'gas.tar.png'
set style fill solid
set multiplot layout 1,1 title "Ions in GO channel"

set yrange [30:90]
set xrange [-60:60]
set xlabel 'X [{\305}]'
set ytics  20
set xtics  20
set ylabel 'Y [{\305}]'
#set size square
set size ratio 0.5
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
p 'gas.tar.dat' u (\$2):(\$3) w lp pt 0 lw 2 lc 3 lt 1  title '' ,\
'OX' u (\$2-${X}):(\$3) w p pt 7 ps 1  lc 2 lt 1  title '',\
'OX' u (\$2):(\$3-${Y}) w p pt 7 ps 1  lc 2 lt 1  title '',\
'OX' u (\$2):(\$3) w p pt 7 ps 1  lc 2 lt 1  title '',\
'OX' u (\$2+${X}):(\$3-${Y}) w p pt 7 ps 1  lc 2 lt 1  title '',\


#'Catom' u (\$2-${X}):(\$3) w p pt 6 ps 2  lc 1 lt 1  title '',\
#'Catom' u (\$2):(\$3-${Y}) w p pt 6 ps 2  lc 1 lt 1  title '',\
#'OX' u (\$2):(\$3) w p pt 7 ps 2  lc 2 lt 1  title '',\
#'OX' u (\$2-${X}):(\$3) w p pt 7 ps 2  lc 2 lt 1  title '',\
#'OX' u (\$2):(\$3-${Y}) w p pt 7 ps 2  lc 2 lt 1  title ''



### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt

