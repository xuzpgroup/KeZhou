#!/bin/bash
shifts=0.89
shifts=0.00
path_G="/home/zhou/WORK_FILE/18-09-02-Ions-EMF/5.Epoxy.new/0.Na.graphene/3.K.benckmark/1.Eq.K.pv-kmesh-9/1.Gr.DOS/LORBIT-11/IPDOS_C.dat"
## plot figure 
cat << EOF > plot.plt
set terminal postscript eps enhanced color  font 'Arial,32' linewidth 1 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'DOS.Int.C.Na.eps'
set style fill solid
set multiplot layout 1,1 title ""
set size ratio 0.72

load "mycol.pal"

# set lmargin 8 
# set rmargin 2 
set key font ',25'
# set key b
set key spacing 1.2

set xlabel '{/Arial-Italic E} (eV)'
# set ytics  2
# set xtics  4
# set xrange [2.0:8.0]
# set yrange [0:40]
set yrange [:40]
set xrange [-30:5]
set ytics  10
set xtics  10
set mxtics  2
set ylabel 'Intg(DOS) (a.u.)'

set key l
# p 'DOS0' u (\$1*1):(\$2)  w l ls 21 lw 4  title '',\



p '< paste IPDOS_C_UP.dat IPDOS_C_DW.dat' \
u (\$1*1+${shifts}):(\$2-\$13)  \
w l ls 21 lw 4  title 'Na|G, s',\
'< paste IPDOS_C_UP.dat IPDOS_C_DW.dat' \
u (\$1*1+${shifts}):(\$5-\$16) \
w l ls 12 lw 4  title '{p_x}',\
'< paste IPDOS_C_UP.dat IPDOS_C_DW.dat' \
u (\$1*1+${shifts}):(\$3-\$14) \
w l ls 13 lw 4  title '{p_y}',\
'< paste IPDOS_C_UP.dat IPDOS_C_DW.dat' \
u (\$1*1+${shifts}):(\$4-\$15) \
w l ls 14 lw 4  title '{p_z}',\
'${path_G}' u (\$1*1):(\$2)  w l ls 21 lw 4 dt 2 title '',\
'${path_G}' u (\$1*1):(\$5)  w l ls 12 lw 4 dt 2 title '',\
'${path_G}' u (\$1*1):(\$3)  w l ls 13 lw 4 dt 2 title '',\
'${path_G}' u (\$1*1):(\$4)  w l ls 14 lw 4 dt 2 title '',\




EOF
gnuplot plot.plt
rm plot.plt

file="DOS.vaspkit"
# convert -density 500 ${file}.eps G.DOS.orbital.png 








