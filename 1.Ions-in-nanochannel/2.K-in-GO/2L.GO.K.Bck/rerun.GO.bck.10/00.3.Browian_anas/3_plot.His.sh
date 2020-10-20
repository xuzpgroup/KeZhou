cat << EOF > plot.plt
set term png  enhanced  font "arial,15" fontscale 1.0  ##size 512, 280 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'His.Lt.png'
set style fill solid
set multiplot layout 1,1 title "interval 10 ps"

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
p 'His.Lt.10ps.dat' u (\$1):(\$2) w p pt 7 ps 1 lw 1 lc 1 lt 1  title 'dx',\
1.01061*exp(-x*x/1.13167) w l lt 1 lc 1 lw 2  title 'Gauss fit',\
'His.Lt.10ps.dat' u (\$3):(\$4) w p pt 7 ps 1 lw 1 lc 2 lt 1  title 'dy',\
0.975333*exp(-x*x/1.21599) w l lt 1 lc 2 lw 2  title 'Gauss fit'



EOF
gnuplot plot.plt
rm plot.plt


