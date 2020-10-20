rm msd.cal
gcc -o msd.cal msd.cal.c  -lm 
./msd.cal

### plot 
#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'VACF.ave.cases.png'
set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set yrange [0:3.2]
set xlabel 't [ps]'
#set ytics  5
set xtics  50
set xrange [0:100]
set ylabel 'D [10^-^9m^2s^-^1]'
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position

set style fill transparent solid 0.2 noborder

p 'msd.ave.dat' u (\$1):(\$4):(\$5) w filledcurves title 'std error' ,\
'msd.ave.dat' u (\$1):(\$2) w lp pt 7 lw 4 lc 1 lt 1  title 'ave',\




 

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt

cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'VACF.all.cases.png'
set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set yrange [0:3.2]
set xlabel 't [ps]'
#set ytics  5
set xtics  50
set xrange [0:100]
set ylabel 'D [10^-^9m^2s^-^1]'
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
#

p  'msd.all.data' u 1:2 w lp lt 1 pt 0 lw 2 lc 1  title 'K^+, case 1' ,\
'msd.all.data' u 1:3 w lp lt 1 pt 0 lw 2 lc 2  title 'case 2' ,\
'msd.all.data' u 1:4 w lp lt 1 pt 0 lw 2 lc 3  title 'case 3' ,\
'msd.all.data' u 1:5 w lp lt 1 pt 0 lw 2 lc 4  title 'case 4' ,\
'msd.all.data' u 1:6 w lp lt 1 pt 0 lw 2 lc 5  title 'case 5' ,\
'msd.all.data' u 1:7 w lp lt 1 pt 0 lw 2 lc 6  title 'case 6' ,\
'msd.all.data' u 1:8 w lp lt 1 pt 0 lw 2 lc 7  title 'case 7' ,\
'msd.all.data' u 1:9 w lp lt 1 pt 0 lw 2 lc 8  title 'case 8' ,\
'msd.all.data' u 1:10 w lp lt 1 pt 0 lw 2 lc 9  title 'case 9' ,\
'msd.ave.dat' u (\$1):(\$2) w lp pt 0 lw 4 lc 0 lt 1  title 'ave',\

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt
