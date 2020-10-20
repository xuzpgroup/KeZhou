#!/bin/bash
## plot figure 
KJ2eV=0.0103643
MidP=24.75
KJmol2kbT=2.493
cat << EOF > plot.plt
# set term post enhanced 20 color eps
# set encoding iso_8859_1
# set term png  enhanced
set terminal postscript eps enhanced color  font 'Arial,24' linewidth 1 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'FES.H10.M1.Conv.eps'
set style fill solid
set multiplot layout 1,1 title ""


#set xlabel '{/Times-Italic d} [{\305}]'
set xlabel '{/Arial-Italic x} [{\305}]'
# set ytics  0.5
# set xtics  5
set xrange [-10:10.0]
set key b
# set yrange [-0.5:4]
set ytics  0.5
# set xtics  1
#set ylabel '{/Times-Italic F} [kJ/mol]'
#set ylabel '{/Times-Italic F} [eV]'
# set ylabel 'F [KJ/mol]'
set ylabel 'PMF [ {/Arial-Italic k}_B{/Arial-Italic T} ]'
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
#set logscale y 
# unset label
# set label 1 at 7, 5.2 'Bulk SPC/E, our = 2.75'
# set label 2 at 7, 5.5 'Bulk TIP3P, PCCP = 6.22'
# set arrow  from 6.0, 2.75 to 24.0, 2.75 nohead lt 2 lc rgb '#DC143C'  ## ad a line

# set style fill transparent solid 0.4 noborder
# set y2range [0:8.0]
# set x2range [-1:1.8]
# set y2label 'N_c'
# set y2tics 2
# set ytics nomirror

# smooth 
# p 'free.1.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 1 lt 1  title '1 ns' smooth  sbezier,\
# 'free.2.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 2 lt 1  title '2 ns' smooth  sbezier,\

p 'free.3.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 3 lt 1  title '3 ns' smooth  sbezier,\
'free.4.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 4 lt 1  title '4 ns' smooth  sbezier,\
'free.5.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 5 lt 1  title '5 ns' smooth  sbezier,\
'free.6.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 6 lt 1  title '6 ns' smooth  sbezier,\
'free.7.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 7 lt 1  title '7 ns' smooth  sbezier,\
'free.8.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 8 lt 1  title '8 ns' smooth  sbezier,\
'free.9.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 9 lt 1  title '9 ns' smooth  sbezier,\
'free.10.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 10 lt 1  title '10 ns' smooth  sbezier,\

# # # raw
# p 'free.1.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 1 lt 1  title '1 ns',\
# 'free.2.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 2 lt 1  title '2 ns' ,\
# 'free.3.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 3 lt 1  title '3 ns' ,\
# 'free.4.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 4 lt 1  title '4 ns' ,\
# 'free.5.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 5 lt 1  title '5 ns' ,\
# 'free.6.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 0 lw 8 lc 6 lt 1  title '6 ns',\



# 'free.out' u (\$1*10):(\$2/${KJmol2kbT}) w lp pt 7 lw 8 lc 2 lt 1  title '' smooth  sbezier


# p 'log' u (\$1*10):(\$2) w l lw 1 lc 1 lt 1  title ''

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt










