rm Cap.cal
gcc -o Cap.cal Cap.cal.c  -lm 
./Cap.cal

### plot 
#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Capture.ave.png'
set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set yrange [0:3.2]
set xlabel 't [ps]'
#set ytics  5
set xtics  250
set yrange [-0.2:1]
set ylabel '<v(0)v(t)>/<v(0)v(0)> '
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position

set style fill transparent solid 0.2 noborder

p 'Cap.ave.dat' u (\$1):(\$4):(\$5) w filledcurves title 'std error' ,\
'Cap.ave.dat' u (\$1):(\$2) w lp pt 0 lw 2 lc 1 lt 1  title 'ave',\




 

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt

cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Cap.all.cases.png'
set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set yrange [0:3.2]
set xlabel 't [ps]'
#set ytics  5
# set xtics  250
# set xrange [0:1000]
set ylabel '<v(0)v(t)>/<v(0)v(0)> '
#set ytics  0.4
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
#

p  '01_norm_rtcf.all.data' u 1:2 w lp lt 1 pt 0 lw 2 lc 1  title 'Na^+, case 1' ,\
'01_norm_rtcf.all.data' u 1:3 w lp lt 1 pt 0 lw 2 lc 2  title 'case 2' ,\
'01_norm_rtcf.all.data' u 1:4 w lp lt 1 pt 0 lw 2 lc 3  title 'case 3' ,\
'01_norm_rtcf.all.data' u 1:5 w lp lt 1 pt 0 lw 2 lc 4  title 'case 4' ,\
'01_norm_rtcf.all.data' u 1:6 w lp lt 1 pt 0 lw 2 lc 5  title 'case 5' ,\
'01_norm_rtcf.all.data' u 1:7 w lp lt 1 pt 0 lw 2 lc 6  title 'case 6' ,\
'01_norm_rtcf.all.data' u 1:8 w lp lt 1 pt 0 lw 2 lc 7  title 'case 7' ,\
'01_norm_rtcf.all.data' u 1:9 w lp lt 1 pt 0 lw 2 lc 8  title 'case 8' ,\
'01_norm_rtcf.all.data' u 1:10 w lp lt 1 pt 0 lw 2 lc 9  title 'case 9' ,\
'Cap.ave.dat' u (\$1):(\$2) w lp pt 0 lw 4 lc 0 lt 1  title 'ave',\

### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt
