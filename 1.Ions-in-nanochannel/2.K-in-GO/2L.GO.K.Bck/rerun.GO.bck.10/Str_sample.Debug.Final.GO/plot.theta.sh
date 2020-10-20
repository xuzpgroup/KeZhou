
cat << EOF > plot.plt
#set term png  enhanced 
set terminal pngcairo  enhanced font "arial,15" fontscale 1.0  ##size 512, 280 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Theta.12.dis.png'
set style fill solid
set multiplot layout 2,1 title ""

# set lmargin 8
# set rmargin 8
# set tmargin 8
# set bmargin 3


set ylabel 'p(cos({/Symbol q})) [a.u.]'
set ylabel 'p({/Symbol q}) [a.u.]'
set xlabel '{/Symbol q}'
#unset ytics
set xrange [0.0:180.0]
#set yrange [0:0.04]
set xtics  30

set arrow  from 0.0, 0.005555 to 180.0, 0.005555 nohead lt 3 lc 0 lw 2  ## ad a line

p 'his.theta.dat' u (\$1):(\$2) w l  lw 3 lc 1 lt 1  title '{/Symbol q}_1',\
'his.theta.dat' u (\$1):(\$3) w l  lw 3 lc 2 lt 1  title '{/Symbol q}_2',\


EOF
gnuplot plot.plt
rm plot.plt


cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Theta.z.dis.png'
set style fill solid
set multiplot layout 2,1 title ""

# set lmargin 8
# set rmargin 8
# set tmargin 8
# set bmargin 3


set ylabel 'p(cos({/Symbol q})) [a.u.]'
set ylabel 'p({/Symbol q}) [a.u.]'
set xlabel '{/Symbol q}'
#unset ytics
set xrange [0.0:90.0]
#set yrange [0:0.04]
set xtics  30
p 'his.thetaz.dat' u (\$1):(\$2) w l  lw 3 lc 1 lt 1  title '{/Symbol q}_z',\


EOF
gnuplot plot.plt
rm plot.plt
