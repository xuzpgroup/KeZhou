cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'His.Lt.png'
set style fill solid
set multiplot layout 1,1 title ""

# set lmargin 8
# set rmargin 8
# set tmargin 8
# set bmargin 3


set ylabel 'Distribution'
set xlabel 'Diffusion length [{\305}]'
#unset ytics
#set xrange [0.0:180.0]
#set yrange [0:10]
#set ytics  2
p 'His.Lt.20ps.dat' u (\$1):(\$2) w p pt 7 ps 1 lw 1 lc 1 lt 1  title 'dx',\
2.27481*exp(-x*x/0.220607) w l lt 1 lc 1 lw 2  title 'Gauss fit',\
'His.Lt.20ps.dat' u (\$3):(\$4) w p pt 7 ps 1 lw 1 lc 2 lt 1  title 'dy',\
2.36568*exp(-x*x/0.201911) w l lt 1 lc 2 lw 2  title 'Gauss fit'

EOF
gnuplot plot.plt
rm plot.plt


