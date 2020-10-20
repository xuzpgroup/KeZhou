#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'pmf.png'
set style fill solid
set multiplot layout 2,1 title ""

#set format x "%.e "
# set xrange [0:12]
# set yrange [0:12]
set xlabel 'r [{\305}]'
# set ytics  2
# set xtics  2
set ylabel 'PMF (kbT)'
#set ytics  0.4
#unset xlabel
#unset xtics
set grid
#set key center at 0.55, 0.8  ## set legend position
p 'rdf.self.dat' u (\$1):(\$4) w lp pt 0 lw 2 lc 2 lt 1  title '',\


#p 'rdf.self.dat' u (\$1):(\$3) w lp pt 0 lw 2 lc 1 lt 1  title 'coord',\



#'rdf_H2O_#2_H2O_[#2r_#2o].csv' u (\$1*1.0E-3):(\$3) w lp pt 0 lw 2 lc 3 lt 1  title 'coord' smooth  sbezier,\


EOF
gnuplot plot.plt
rm plot.plt


