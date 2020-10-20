
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Pair.t.png'
set style fill solid
set multiplot layout 1,1 title ""

# set lmargin 8
# set rmargin 8
# set tmargin 8
# set bmargin 3


set ylabel 'N_{pair}'
set xlabel 't [ps]'
#unset ytics
set xrange [0.0:10.0]
# set yrange [0:10]
# set ytics  2

p 'pair.t.dat' u (\$1/1000):(\$3) w l  lw 3 lc 1 lt 1  title 'N_c'

EOF
gnuplot plot.plt
rm plot.plt


