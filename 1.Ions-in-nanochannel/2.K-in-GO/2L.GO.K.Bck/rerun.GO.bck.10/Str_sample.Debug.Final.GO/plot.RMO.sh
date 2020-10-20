
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'His.Rmo.png'
set style fill solid
set multiplot layout 1,1 title ""

# set lmargin 8
# set rmargin 8
# set tmargin 8
# set bmargin 3


set ylabel 'P(r)'
set xlabel 'r [{\305}]'
#unset ytics
#set xrange [0.0:180.0]
set xrange [0:3.5]
set xtics  0.5
p 'his.Rmo.dat' u (\$1):(\$2) w lp  lw 3 pt 4 lc 1 lt 1  title 'R_M_O',\
'his.Rmo.dat' u (\$1):(\$3) w lp  lw 3 pt 4 lc 2 lt 1  title 'R_O_H',\


EOF
gnuplot plot.plt
rm plot.plt


