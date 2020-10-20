#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Nc.his.png'
set style fill solid
set multiplot layout 1,1 title ""


set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
#set xtic rotate by -45 scale 0

#set format x "%.e "
#set xtics 250
#set ytics 10
#set xtics  2
#set xrange [1:8]
#set yrange [3:-12]

set ylabel 'P(N_c) [%]'
#set xlabel 'Snopshot'
#set xtics  5
#set ytics  0.005
#set ytics 0.1

#set ytics 0.1
p 'his.nc.dat' u 2:xtic(1)  title ''  #ti col 

#p '1_CO2/Bind.dat' u (\$1*1.0E-5-2):(\$4-\$7-\$10) w p pt 5 ps 2 lw 2 lc 1 lt 1  title 'CO_2',\
#'2_N2/Bind.dat' u (\$1*1.0E-5-2):(\$4-\$7-\$10) w p pt 7 ps 2 lw 2 lc 2 lt 1  title 'N_2',\
#'3_H2/Bind.dat' u (\$1*1.0E-5-2):(\$4-\$7-\$10) w p pt 9 ps 2 lw 2 lc 3 lt 1  title 'H_2',\
#'4_CH4/Bind.dat' u (\$1*1.0E-5-2):(\$4-\$7-\$10) w p pt 11 ps 2 lw 2 lc 4 lt 1  title 'CH_4',\



EOF
gnuplot plot.plt
rm plot.plt









