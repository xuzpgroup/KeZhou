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
set output 'ACF.single.png'
set style fill solid
set multiplot layout 1,1 title ""


set xrange [0:10]
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
p './VACF.single/VACF.single.1.data' u (\$1):(\$2+\$3) \
w l lw 2 lc 1 lt 1  title '1',\
'./VACF.single/VACF.single.2.data' u (\$1):(\$2+\$3) \
w l lw 2 lc 2 lt 1  title '2',\
'./VACF.single/VACF.single.3.data' u (\$1):(\$2+\$3) \
w l lw 2 lc 3 lt 1  title '3',\
'./VACF.single/VACF.single.3.data' u (\$1):(\$2+\$3) \
w l lw 2 lc 4 lt 1  title '4',\
'./VACF.single/VACF.single.31.data' u (\$1):(\$2+\$3) \
w l lw 2 lc 1 lt 1  title '31',\
'./VACF.single/VACF.single.32.data' u (\$1):(\$2+\$3) \
w l lw 2 lc 2 lt 1  title '32',\
'./VACF.single/VACF.single.3.data' u (\$1):(\$2+\$3) \
w l lw 2 lc 3 lt 1  title '33',\
'./VACF.single/VACF.single.3.data' u (\$1):(\$2+\$3) \
w l lw 2 lc 4 lt 1  title '34',\


EOF
gnuplot plot.plt
rm plot.plt

