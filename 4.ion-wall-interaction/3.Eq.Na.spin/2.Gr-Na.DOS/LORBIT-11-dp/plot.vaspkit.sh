#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set terminal postscript eps enhanced color  font 'Arial,32' linewidth 1 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'DOS.vaspkit.K.eps'
set style fill solid
set multiplot layout 1,1 title ""
set size ratio 1.4

load "mycol.pal"

# set lmargin 8 
# set rmargin 2 
set key font ',30'
# set key b

set xlabel '{/Arial-Italic E} (eV)'
# set ytics  2
# set xtics  4
set xrange [-5:2.0]
set xtics  2
set mxtics  2
set yrange [-20:20]
set ytics  10
# set xtics  10
# set ylabel 'DOS (states/eV/vell)'
set ylabel 'DOS (a.u.)'
unset xlabel

set key r
set key spacing 1.2

set arrow  from 0.0, -20 to 0.0, 20 nohead dt 2 lw 2 lc rgb 'blue'  ## ad a line


# p 'TDOS.dat' u (\$1*1):(\$2)  w l ls 21 lw 4  title 'DOS',\
# 'TDOS.dat' u (\$1*1):(\$3)  w l ls 21 lw 4  title '',\
# 'PDOS_Cu_UP.dat' u (\$1*1):(\$2)  w l ls 22 lw 4  title 'Cu, s',\
# 'PDOS_Cu_DW.dat' u (\$1*1):(\$2)  w l ls 22 lw 4  title '',\
# 'PDOS_Cu_UP.dat' u (\$1*1):(\$3+\$4+\$5)  w l ls 23 lw 4  title 'Cu, p',\
# 'PDOS_Cu_DW.dat' u (\$1*1):(\$3+\$4+\$5)  w l ls 23 lw 4  title '',\
# # 'PDOS_Cu_UP.dat' u (\$1*1):(\$6+\$7+\$8+\$9+\$10)  w l ls 24 lw 4  title 'Cu, d',\
# # 'PDOS_Cu_DW.dat' u (\$1*1):(\$6+\$7+\$8+\$9+\$10)  w l ls 24 lw 4  title '',\


p 'TDOS.dat' u (\$1*1):(\$2)  w l ls 21 lw 4  title 'DOS',\
'TDOS.dat' u (\$1*1):(\$3)  w l ls 21 lw 4  title '',\
'PDOS_C_UP.dat' u (\$1*1):(\$4)  w l ls 23 lw 4  title '{/Arial-Italic p_z}',\
'PDOS_C_DW.dat' u (\$1*1):(\$4)  w l ls 23 lw 4  title '',\
'PDOS_Na_UP.dat' u (\$1*1):(\$2*2)  w l ls 22 lw 4  title '{/Arial-Italic s}',\
'PDOS_Na_DW.dat' u (\$1*1):(\$2*2)  w l ls 22 lw 4  title '',\




# 'PDOS_Cu_UP.dat' u (\$1*1):(\$11)  w l ls 22 lw 4  title 'Cu',\
# 'PDOS_Cu_DW.dat' u (\$1*1):(\$11)  w l ls 22 lw 4  title '',\

# 'PDOS_C_UP.dat' u (\$1*1):(\$11)  w l ls 13 lw 2  title 'G, up',\
# 'PDOS_C_DW.dat' u (\$1*1):(\$11)  w l ls 14 lw 2  title 'G, down',\



EOF
gnuplot plot.plt
rm plot.plt

file="DOS.vaspkit"
convert -density 500 ${file}.eps G-Cu.DOS.orbital.png 








