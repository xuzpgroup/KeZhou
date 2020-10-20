E0=-9.22077616
cat << EOF > plot.plt
# set term post enhanced 20 color eps
# set encoding iso_8859_1
set terminal postscript eps enhanced color  font 'Arial,25' linewidth 1 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'WorkFunction.eps'
set style fill solid
set multiplot layout 1,1 title ""
# set size ratio 1.5
set size ratio 0.72



#set xlabel '{/Times-Italic d} [{\305}]'
set xlabel '{/Arial-Italic z} ({\305})'
# set ytics  0.1
set xtics  5
# set format x "%.1f"
# set format y "%.1f"
# set xrange [6.0:24.0]
# set yrange [:0]
#set ytics  5
# set xtics  1
#set ylabel '{/Times-Italic F} [kJ/mol]'
#set ylabel '{/Times-Italic F} [eV]'
set ylabel '{/Symbol F} (eV)'
# set ylabel 'PMF [ {/Arial-Italic k}_B{/Arial-Italic T} ]'
set ytics  5
#unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
#set logscale y 
# unset label

# set arrow  from -10.0, 0.0 to  10.0, 0.0 nohead dt 1 lw 2 lc rgb 'black'  ###rgb '#DC143C' 
# set arrow  from  2.2, -0.4 to   2.2, 0.8 nohead dt 2 lc rgb 'red'  ###rgb '#DC143C' 


p 'vplanar.txt' u (\$2*15):(\$3) w l lw 4 lc 0 lt 1  title '',\
'vplanar.txt' u ((\$2+1.0)*15):(\$3) w l lw 4 lc 0 lt 1  title '',\
'vplanar.txt' u ((\$2-1.0)*15):(\$3) w l lw 4 lc 0 lt 1  title '',\


### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt

# convert -density 500 Ed.eps Ed.png 

