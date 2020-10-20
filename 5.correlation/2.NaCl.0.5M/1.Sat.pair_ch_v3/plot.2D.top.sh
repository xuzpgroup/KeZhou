#!/bin/bash
## plot figure 
cat << EOF > plot.plt
# set term png  enhanced 
set terminal postscript eps enhanced color  font 'Arial,25' linewidth 1
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Top.2D.eps'

set style fill solid
set multiplot layout 1,1 title ""
set size ratio 1.0

# set ylabel '{/Arial-Italic y} ({\305})'
# set xlabel '{/Arial-Italic x} ({\305})'
set xlabel '{/Arial-Italic n}_{+}'
set ylabel '{/Arial-Italic n}_{-}'
#set ytics 0.1
# load 'Moreland.plt'
# load 'plasma.pal'
# load 'parula.pal'
# load 'viridis.pal'

# unset ylabel
# unset ytics
set cbrange [-2:2]
# set cblabel "{/Symbol D}{/Arial-Italic G} ( {/Arial-Italic k}_B{/Arial-Italic T} )"
# set cblabel "log {/Arial-Italic n}_{ij}"
set cblabel "{/Arial-Italic n}_{ij}"

set yrange [0:11]
set xrange [0:11]

# set palette defined (0 "white", 1 "blue", 2 "cyan", 3 "red", 4 "orange")
set palette defined (0 "white", 1 "blue", 2 "red")
# set palette defined (0 "green", 1 "blue", 2 "cyan", 3 "red", 4 "orange")

set xtics 2 offset 0.5,0
set ytics 2 offset 0,0.5
set mxtics 2
set mytics 2
set pm3d map
set view map scale 1
set grid xtics ytics mxtics mytics ls 1 lc 0
# set xtics offset 0.5
# set ytics offset 0.5
# unset xtics
# set xtics ("0" 0.5, "2" 2.5, "4" 4.5, "6" 6.5, "8" 8.5, "10" 10.5)
set arrow  from 0.0, 0.0 to 11.0, 11.0 nohead lt 1 lc 0  ## ad a line

set cbtics ("10^{-2}" -2, "10^{-1}" -1, "1" 0, "10^1" 1, "10^2" 2)


sp 'cluster.top.dat' u (\$1+0.5):(\$2+0.5):(log(\$3+0.0001)/log(10.0))\
  w image ti '',\





EOF
gnuplot plot.plt
rm plot.plt


# convert -density 500 Pot.2D.his.eps Pot.2D.his.png 
