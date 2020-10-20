E0=-9.22077616
cat << EOF > plot.plt
# set term post enhanced 20 color eps
# set encoding iso_8859_1
set terminal postscript eps enhanced color  font 'Arial,25' linewidth 1 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'dch.density.Na.eps'
set style fill solid
set multiplot layout 1,1 title ""
# set size ratio 1.5
set size ratio 0.3



#set xlabel '{/Times-Italic d} [{\305}]'
set xlabel '{/Arial-Italic z} ({\305})'
# set ytics  0.1
set xtics  2
# set format x "%.1f"
# set format y "%.1f"
set xrange [-8:8]
set yrange [-0.6:0.6]
set ytics  0.3
set xtics  4
#set ylabel '{/Times-Italic F} [kJ/mol]'
#set ylabel '{/Times-Italic F} [eV]'
set ylabel '{/Symbol D}{/Symbol-Italic r}({/Arial-Italic z}) (e/{\305})'
# set ylabel 'PMF [ {/Arial-Italic k}_B{/Arial-Italic T} ]'
#set ytics  0.4
unset xlabel
#unset xtics
#set key center at 0.55, 0.8  ## set legend position
#set logscale y 
# unset label

set arrow  from -8.0,  0.0 to    8.0, 0.0 nohead dt 1 lw 2 lc rgb 'black'  ###rgb '#DC143C' 
set arrow  from  2.30, -0.6 to   2.30, 0.6 nohead dt 2 lc rgb 'red'  ###rgb '#DC143C' 


p './His.in.z.dat' u (\$1-10):(\$2) w l lw 4 lc 0 lt 1  title 'Na',\


### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt

convert -density 500 dch.density.Na.eps dch.density.Na.png 

