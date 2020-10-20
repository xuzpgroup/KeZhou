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
set output 'RTCF.png'
set style fill solid
set multiplot layout 1,1 title ""


#set xrange [0:1000]
#set yrange [-2:6]
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
# p '0rtcf.data' u (\$1):(\$2) w lp pt 2 lw 2 lc 1 lt 1  title ''

p '01_norm_rtcf.data' u (\$1):(\$2) w lp pt 6 lw 2 lc 1 lt 1  title '',\
 0.631038*exp(-x/243.613)+(1-0.631038)*exp(-x/20.1103) w l lt 1 lc 2 lw 2  title ""

# p '0rtcf.data' u (\$1):(\$2) w l lw 2 lc 1 lt 1  title '',\

# exp(-x/27.5513) w l lt 1 lc 2 lw 2  title "",\
# exp(-x/16.407) w l lt 1 lc 2 lw 2  title "",\
# exp(-x/1.407) w l lt 1 lc 2 lw 2  title "",\



EOF
gnuplot plot.plt
rm plot.plt


