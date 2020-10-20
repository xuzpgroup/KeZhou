Cal2J=4.184
#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1
#set title "Time ${t} ps"
set output 'Drag.ion.PMF.png'
set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set xtics 250
#set ytics 10
#set xtics  0.5
#set xrange [0.3:3.2]
#set yrange [0:60]
#set xlabel 'Time/ns'
#set xtics  0.3
#set ytics  20
set xlabel 'X_I_o_n [{\305}]'
set ylabel 'PMF [kj/mol]'
#set xrange [0.1:0.5]
#set ytics  200
p 'smd.pmf.data' u (\$1*1.0E-6*0.25):(\$8*${Cal2J}) w lp pt 0 lw 2 lc 1 lt 1  title ''


EOF
gnuplot plot.plt
rm plot.plt

