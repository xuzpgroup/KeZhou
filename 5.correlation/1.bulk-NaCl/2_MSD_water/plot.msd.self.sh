#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'msd.self.png'
set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set xrange [0:5.0]
set xlabel 't [ps]'
#set ytics  5
#set xtics  100
set ylabel 'MSD [nm^2]'
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
p 'msd.data' u (\$1):(\$2/100) w lp pt 4 lw 2 lc 3 lt 1  title '' 

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt


#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'diffusion.self.png'
set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set xtics 250
#set ytics 10
#set xrange [0:850]
#set yrange [1.0:1.4]
#set ylabel 'Density (gcm^-^3)'
#set xlabel 'Time/ps'
#set xtics  5
#set ytics  0.005
#set ytics 0.1

#set ytics 0.1
#set xrange [0.2:1.0]
#set xrange [0.0:1.0]
#set yrange [1:3.9]
#set yrange [0:200]
#set xrange [0:250]
#set yrange [0:4.5]
set xlabel 't [ps]'
#set ytics  0.5
#set xtics  50
set ylabel 'D [10^-^9m^2s^-^1]'
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
p 'msd.data' u (\$1):(\$3) w lp pt 0 lw 2 lc 3 lt 1  title '' 

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt


