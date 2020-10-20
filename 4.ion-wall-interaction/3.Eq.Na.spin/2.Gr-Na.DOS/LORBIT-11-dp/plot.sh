#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set terminal postscript eps enhanced color  font 'Arial,25' linewidth 1 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'DOS.eps'
set style fill solid
set multiplot layout 1,1 title ""
set size ratio 0.72

load "mycol.pal"

# set lmargin 8 
# set rmargin 2 
set key font ',20'
# set key b

set xlabel '{/Arial-Italic E} (eV)'
# set ytics  2
# set xtics  4
set xrange [-25:10.0]
set xtics  5
# set yrange [0:40]
set ytics  10
# set xtics  10
# set ylabel 'DOS (states/eV/vell)'
set ylabel 'DOS (a.u.)'

set key l

set arrow  from 0.0, -20 to 0.0, 20 nohead dt 2 lw 2 lc rgb '#DC143C'  ## ad a line


p 'TDOS.dat' u (\$1*1):(\$2)  w l ls 21 lw 4  title 'Na/G, DOS',\
'TDOS.dat' u (\$1*1):(\$3)  w l ls 21 lw 4  title '',\
'PDOS_C_UP.dat' u (\$1*1):(\$2)  w l ls 11 lw 4  title 'C, {/Arial-Italic s}',\
'PDOS_C_DW.dat' u (\$1*1):(\$2)  w l ls 11 lw 4  title '',\
'PDOS_C_UP.dat' u (\$1*1):(\$5)  w l ls 12 lw 4  title 'C, {/Arial-Italic p_x}',\
'PDOS_C_DW.dat' u (\$1*1):(\$5)  w l ls 12 lw 4  title '',\
'PDOS_C_UP.dat' u (\$1*1):(\$3)  w l ls 13 lw 4  title 'C, {/Arial-Italic p_y}',\
'PDOS_C_DW.dat' u (\$1*1):(\$3)  w l ls 13 lw 4  title '',\
'PDOS_C_UP.dat' u (\$1*1):(\$4)  w l ls 14 lw 4  title 'C, {/Arial-Italic p_z}',\
'PDOS_C_DW.dat' u (\$1*1):(\$4)  w l ls 14 lw 4  title '',\

# 'PDOS_C_UP.dat' u (\$1*1):(\$6+\$7+\$8+\$9+\$10)  w l ls 13 lw 4  title 'C, {/Arial-Italic d}',\
# 'PDOS_C_DW.dat' u (\$1*1):(\$6+\$7+\$8+\$9+\$10)  w l ls 13 lw 4  title '',\





EOF
gnuplot plot.plt
rm plot.plt

file="DOS"
convert -density 500 ${file}.eps ${file}.Na.png 








