
#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1
#set title "Time ${t} ps"
set output 'Water.z.O.png'
set style fill solid
set multiplot layout 2,1 title ""

set xlabel 'Position Z({\305})'
#set yrange [0:0.03]
#set xrange [0.0:8.386]
#set xtics 2
#set ytics  0.1
#unset ytics
set style fill solid 0.4

set ylabel 'O atom distribution probability'
#set ytics  200
p 'Water.distribution.data' u (\$1):(\$2) w filledcurves lw 2 lc 7 lt 1  title 'h=7.0{\305}'

#set ylabel 'Mass distribution probability'
#p 'IL.distribution.data' u (\$1-5.0):(\$5) w filledcurves lw 2 lc 7 lt 1  title 'h=11.5{\305}'

EOF
gnuplot plot.plt
rm plot.plt

cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1
#set title "Time ${t} ps"
set output 'Water.z.mass.png'
set style fill solid
set multiplot layout 2,1 title ""

set xlabel 'Position Z({\305})'
#set yrange [0:0.03]
#set xrange [0.0:8.386]
#set xtics 2
#set ytics  0.1
#unset ytics
set style fill solid 0.4

set ylabel 'Mass distribution probability'
#set ytics  200
p 'Water.distribution.data' u (\$1):(\$3) w filledcurves lw 2 lc 7 lt 1  title 'h=7.0{\305}'

#set ylabel 'Mass distribution probability'
#p 'IL.distribution.data' u (\$1-5.0):(\$5) w filledcurves lw 2 lc 7 lt 1  title 'h=11.5{\305}'

EOF
gnuplot plot.plt
rm plot.plt

cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1
#set title "Time ${t} ps"
set output 'Water.z.particle.png'
set style fill solid
set multiplot layout 2,1 title ""

set xlabel 'Position Z({\305})'
#set yrange [0:0.03]
#set xrange [0.0:8.386]
#set xtics 2
#set ytics  0.1
#unset ytics
set style fill solid 0.4

set ylabel 'Atom distribution probability'
#set ytics  200
p 'Water.distribution.data' u (\$1):(\$3) w filledcurves lw 2 lc 7 lt 1  title 'h=7.0{\305}'

#set ylabel 'Mass distribution probability'
#p 'IL.distribution.data' u (\$1-5.0):(\$5) w filledcurves lw 2 lc 7 lt 1  title 'h=11.5{\305}'

EOF
gnuplot plot.plt
rm plot.plt
