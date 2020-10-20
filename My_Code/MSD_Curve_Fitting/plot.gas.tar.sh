X=63.000000
Y=111.544072

#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'gas.tar.png'
set style fill solid
set multiplot layout 1,1 title ""

set yrange [10:70]
set xrange [10:70]
set xlabel 'X [{\305}]'
set ytics  10
set xtics  10
set ylabel 'Y [{\305}]'
set size square
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
p 'gas.tar.dat' u (\$2):(\$3) w lp pt 0 lw 2 lc 3 lt 1  title 'Na^+, h=6.8{\305}, G' ,\
'Background.dat' u (\$2):(\$3) w p pt 6 ps 0 lc 1 lt 1  title '',\
'Background.dat' u (\$2+${X}):(\$3) w p pt 6 ps 0 lc 1 lt 1  title '',\

#'Catom' u (\$2-${X}):(\$3) w p pt 6 ps 2  lc 1 lt 1  title '',\
#'Catom' u (\$2):(\$3-${Y}) w p pt 6 ps 2  lc 1 lt 1  title '',\
#'OX' u (\$2):(\$3) w p pt 7 ps 2  lc 2 lt 1  title '',\
#'OX' u (\$2-${X}):(\$3) w p pt 7 ps 2  lc 2 lt 1  title '',\
#'OX' u (\$2):(\$3-${Y}) w p pt 7 ps 2  lc 2 lt 1  title ''



### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt

