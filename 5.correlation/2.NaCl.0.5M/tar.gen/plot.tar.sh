#!/bin/bash
## plot figure 
H=10.0
Lx=8;
Ly=8;
R=2
cat << EOF > plot.plt
# set term png  enhanced 
set terminal postscript eps enhanced color  font 'Arial,25' linewidth 1
set encoding iso_8859_1
set output 'Tar.ions.50.eps'
set style fill solid
set multiplot layout 1,1 title ""
set size ratio  1.000

#set format x "%.e "
#set xtics 250
set ytics 1
set xtics 1
# set xrange [-${R}:${R}]
# set yrange [-${R}:${R}]

set ylabel '{/Arial-Italic y} (nm)'
set xlabel '{/Arial-Italic x} (nm)'

#set ytics 0.1
# load 'Moreland.plt'
# load 'plasma.pal'
# load 'parula.pal'
# load 'viridis.pal'
load 'mycol.pal'

# p 'Tar.dat' u (\$6/10):(\$7/10) every 2::1 w lp pt 1   lc rgb 'black' ti ""


p 'Tar.dat' u (\$3/10):(\$4/10) every 2::1 w l ls 11   ti "",\
'Tar.dat' u (\$6/10):(\$7/10) every 2::1 w l ls 12    ti "",\


EOF
gnuplot plot.plt
rm plot.plt


