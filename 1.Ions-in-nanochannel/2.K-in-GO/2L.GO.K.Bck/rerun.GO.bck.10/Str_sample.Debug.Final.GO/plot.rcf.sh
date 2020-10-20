
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Residence.cf.png'
set style fill solid
set multiplot layout 1,1 title ""

# set lmargin 8
# set rmargin 8
# set tmargin 8
# set bmargin 3


set ylabel 'Residence correlation function'
set xlabel 'dt [ps]'
#unset ytics
#set xrange [0.0:180.0]
#set xrange [0:10]
#set xtics  2
p 'RCF.dat' u (\$1):(\$2) w lp  lw 3 pt 4 lc 1 lt 1  title '',\


EOF
gnuplot plot.plt
rm plot.plt


