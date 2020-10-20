T=0.75
PI=3.1415926

#!/bin/bash
## plot figure 
cat << EOF > plot.plt
#set terminal postscript eps color solid linewidth 2 "Helvetica" 26
# set terminal pngcairo  transparent enhanced font "arial,12" fontscale 1.0  ##size 512, 280 
set term png  enhanced  font "arial,15" fontscale 1.0  ##size 512, 280 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Velocity.ACF.I.png'
set style fill solid
set multiplot layout 1,1 title ""


set xrange [0:1]
#set yrange [-2:6]
#set yrange [0:5]
set xlabel 't [ps]'
#set ylabel 'Q ({/Symbol-Oblique w})'
set ylabel '<v_M(0)v_C(t)> (A^2/ps^2)'
#set xtics  5
#unset xlabel
#unset xtics
#set key center at 130, 4.0  ## set legend position
#set logscale

#p 'Q.omega.data' u (\$1*2*${PI}/${T}*1.0E-3):(\$3) w p pt 0 lw 2 lc 0 lt 1  title ''
p 'vacf_self_MM.data' u (\$1):(\$5*1.0E-4) w lp pt 6 lw 2 lc 1 lt 1  title 'xyz, Na^+, II',\
'vacf_self_CC.data' u (\$1):(\$5*1.0E-4) w lp pt 6 lw 2 lc 2 lt 1  title 'xyz, Cl^-, II',\
'vacf_cross_MC.data' u (\$1):(\$5*1.0E-4) w lp pt 6 lw 2 lc 3 lt 1  title 'xyz, Na^+-Cl^-, IL',\
'vacf_cross_MM.data' u (\$1):(\$5*1.0E-4) w lp pt 6 lw 2 lc 4 lt 1  title 'xyz, Na^+-Na^+, IL',\
'vacf_cross_CC.data' u (\$1):(\$5*1.0E-4) w lp pt 6 lw 2 lc 5 lt 1  title 'xyz, Cl^--Cl^-, IL',\




EOF
gnuplot plot.plt
rm plot.plt


