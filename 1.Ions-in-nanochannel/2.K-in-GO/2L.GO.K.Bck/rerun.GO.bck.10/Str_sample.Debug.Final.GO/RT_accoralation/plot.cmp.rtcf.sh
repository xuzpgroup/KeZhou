T=0.75
PI=3.1415926

#!/bin/bash
## plot figure 
cat << EOF > plot.plt
#set terminal postscript eps color solid linewidth 2 "Helvetica" 26
# set terminal pngcairo  transparent enhanced font "arial,12" fontscale 1.0  ##size 512, 280 
set terminal pngcairo  enhanced font "arial,15" fontscale 1.0  ##size 512, 280 
# set term png  enhanced  font "arial,15" fontscale 1.0  ##size 512, 280 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'RTCF.all.png'
set style fill solid
set multiplot layout 1,1 title "Residence time correlation function"


set xrange [0:100]
set yrange [0:1.4]
set ytics  0.2
#set yrange [0:5]
set xlabel 't [ps]'
#set ylabel 'Q ({/Symbol-Oblique w})'
set ylabel '<v(0)v(t)> '
#set xtics  5
#unset xlabel
#unset xtics
#set key center at 130, 4.0  ## set legend position
#set logscale

#p 'Q.omega.data' u (\$1*2*${PI}/${T}*1.0E-3):(\$3) w p pt 0 lw 2 lc 0 lt 1  title ''
p '01_norm_rtcf.Na.data' u (\$1):(\$2) w l lw 2 lc 1 lt 1  title 'Na^+',\
'01_norm_rtcf.K.data' u (\$1):(\$2) w l lw 2 lc 2 lt 1  title 'K^+',\
'01_norm_rtcf.Mg.data' u (\$1):(\$2) w l lw 2 lc 3 lt 1  title 'Mg^2^+',\
'01_norm_rtcf.Ca.data' u (\$1):(\$2) w l lw 2 lc 4 lt 1  title 'Ca^2^+',\
'01_norm_rtcf.Cl.data' u (\$1):(\$2) w l lw 2 lc 5 lt 1  title 'Cl^-',\


EOF
gnuplot plot.plt
rm plot.plt


