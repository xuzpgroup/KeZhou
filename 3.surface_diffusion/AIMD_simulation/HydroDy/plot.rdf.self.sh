#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'rdf.self.png'
set style fill solid
set multiplot layout 2,1 title ""

#set format x "%.e "
set xrange [0:1]
# set yrange [0:4.0]
set xlabel 'r (nm)'
set ytics 1
set xtics  0.2
set mxtics  2
set ylabel 'g(r)'
#set ytics  0.4
#unset xlabel
#unset xtics
#set grid
#set key center at 0.55, 0.8  ## set legend position
p 'rdf.8.8.dat' u (\$1*.1):(\$2) w lp pt 0 lw 2 lc 1 lt 1  title 'rdf',\

set yrange [0:10]
p 'rdf.8.8.dat' u (\$1*.1):(\$3) w lp pt 0 lw 2 lc 1 lt 1  title 'N_c',\

#'rdf.self.dat' u (\$1*.1):(\$3) w lp pt 0 lw 2 lc 1 lt 1  title 'coord',\



#'rdf_H2O_#2_H2O_[#2r_#2o].csv' u (\$1*1.0E-3):(\$3) w lp pt 0 lw 2 lc 3 lt 1  title 'coord' smooth  sbezier,\


EOF
gnuplot plot.plt
rm plot.plt


