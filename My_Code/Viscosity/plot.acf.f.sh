T=0.75
PI=3.1415926

#!/bin/bash
## plot figure 
cat << EOF > plot.plt
#set terminal postscript eps color solid linewidth 2 "Helvetica" 26
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Press.ACF.png'
set style fill solid
set multiplot layout 1,1 title ""


set xrange [0:2400]
#set yrange [-2:6]
#set yrange [0:5]
set xlabel 't [ps]'
#set ylabel 'Q ({/Symbol-Oblique w})'
set ylabel '<P(0)P(t)> '
#set xtics  5
#unset xlabel
#unset xtics
#set key center at 130, 4.0  ## set legend position
#set logscale

#p 'Q.omega.data' u (\$1*2*${PI}/${T}*1.0E-3):(\$3) w p pt 0 lw 2 lc 0 lt 1  title ''
p '01_norm_facf.data' u (\$1*1000):(\$2) w lp pt 0 lw 2 lc 1 lt 1  title 'p_x_y',\
'01_norm_facf.data' u (\$1*1000):(\$3) w lp pt 0 lw 2 lc 2 lt 1  title 'p_x_z',\
'01_norm_facf.data' u (\$1*1000):(\$4) w lp pt 0 lw 2 lc 3 lt 1  title 'p_y_z',\
'01_norm_facf.data' u (\$1*1000):(\$5) w lp pt 0 lw 2 lc 4 lt 1  title 'ave',\


EOF
gnuplot plot.plt
rm plot.plt


cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Press.int.ACF.png'
set style fill solid
set multiplot layout 2,2 title ""



#set yrange [1:3.9]
#set yrange [0:5]
set xlabel 't [ps]'
#set ylabel 'Q ({/Symbol-Oblique w})'
set ylabel 'Viscosity [ mPa.s]'
#set xtics  5
#unset xlabel
#unset xtics
set xtics 200
set xrange [0:1000]
#set key center at 130, 4.0  ## set legend position
#set logscale

#p 'Q.omega.data' u (\$1*2*${PI}/${T}*1.0E-3):(\$3) w p pt 0 lw 2 lc 0 lt 1  title ''
p '03.trap.acf.dat' u (\$1):(\$3*1000) w lp pt 4 lw 2 lc 1 lt 1  title 'xy',\

p '03.trap.acf.dat' u (\$1):(\$4*1000) w lp pt 4 lw 2 lc 2 lt 1  title 'xz',\

p '03.trap.acf.dat' u (\$1):(\$5*1000) w lp pt 4 lw 2 lc 3 lt 1  title 'yz',\

p '03.trap.acf.dat' u (\$1):(\$2*1000) w lp pt 4 lw 2 lc 4 lt 1  title 'ave',\


EOF
gnuplot plot.plt
rm plot.plt


