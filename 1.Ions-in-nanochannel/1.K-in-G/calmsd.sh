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
set output 'msd.ave.cases.png'
set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set yrange [0:3.2]
set xlabel 't [ps]'
#set ytics  5
set xtics  250
set xrange [0:500]
set ylabel 'MSD [{\305}^2]'
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
p 'msd.ave.dat' u (\$1):(\$2) w lp pt 0 lw 4 lc 1 lt 1  title 'ave',\

 

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt

cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'msd.all.cases.png'
set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set yrange [0:3.2]
set xlabel 't [ps]'
#set ytics  5
set xtics  250
set xrange [0:500]
set ylabel 'MSD [{\305}^2]'
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
#

p  for [i=2:11:1] 'msd.all.data' u 1:i w lp lt 1 pt 0 lw 2 lc (i-1)  title 'casr ${i}' ,\
'msd.ave.dat' u (\$1):(\$2) w lp pt 0 lw 4 lc 0 lt 1  title 'ave',\

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt
