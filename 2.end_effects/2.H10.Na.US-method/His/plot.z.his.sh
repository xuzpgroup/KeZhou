
#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1
#set title "Time ${t} ps"
set output 'dis.us.png'
set style fill solid
set multiplot layout 2,1 title ""

set xlabel 'Position Z({\305})'
#set yrange [0:0.03]
# set xrange [-2:15]
#set xtics 2
#set ytics  0.1
#unset ytics
#set style fill solid 0.4

set ylabel 'Count'
#set ytics  200

p 'distribution.data' u (\$1*10):(\$2) w l  lw 2 lc rgb '#FF0000' lt 1  title '',\



#set ylabel 'Mass distribution probability'
#p 'IL.distribution.data' u (\$1-5.0):(\$5) w filledcurves lw 2 lc 7 lt 1  title 'h=11.5{\305}'

EOF
gnuplot plot.plt
rm plot.plt


