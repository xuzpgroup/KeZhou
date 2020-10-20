
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'SPAN.t.png'
set style fill solid
set multiplot layout 1,1 title ""

# set lmargin 8
# set rmargin 8
# set tmargin 8
# set bmargin 3


set ylabel 'P(dr)'
set xlabel 'dr [{\305}]'
#unset ytics
#set xrange [0.0:180.0]
set xrange [0:10]
set xtics  2
p 'his.span.dat' u (\$1):(\$2) w lp  lw 3 pt 4 lc 1 lt 1  title 'd_X',\
'his.span.dat' u (\$1):(\$3) w lp  lw 3 pt 4 lc 2 lt 1  title 'd_Y',\
'his.span.dat' u (\$1):(\$4) w lp pt 4  lw 3 lc 3 lt 1  title 'd_Z',\


EOF
gnuplot plot.plt
rm plot.plt


