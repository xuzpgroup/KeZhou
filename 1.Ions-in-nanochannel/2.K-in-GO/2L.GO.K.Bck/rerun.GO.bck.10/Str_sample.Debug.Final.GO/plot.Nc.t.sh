
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Nc.t.png'
set style fill solid
set multiplot layout 2,1 title ""

# set lmargin 8
# set rmargin 8
# set tmargin 8
# set bmargin 3


set ylabel 'N_c'
set xlabel 't [ps]'
#unset ytics
#set xrange [0.0:180.0]
set yrange [0:10]
set ytics  2
p 'info.Nc.dat' u (\$1):(\$2) w l  lw 3 lc 1 lt 1  title 'N_c'

EOF
gnuplot plot.plt
rm plot.plt


