atm2mpa=0.101324
h=10.2
fG=`echo "scale=3;($h-4.5)/($h-3.4)"|bc`;
fGO=`echo "scale=3;($h-4.5)/($h-4.8)"|bc`;
#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1
#set title "Time ${t} ps"
set output 'eq.T.png'
set style fill solid
set multiplot layout 1,1 title ""

#set yrange [-600:400]
#set xlabel 'Time/ns'
#set xtics  0.3
#set ytics  20
set ylabel 'Temp'
set xlabel 't'
#set ytics  200
p 'log_p' u (\$1*1.0E-6):(\$3) w lp pt 0 lw 2 lc 2 lt 1  title 'T_w',\


EOF
gnuplot plot.plt
rm plot.plt

#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1
#set title "Time ${t} ps"
set output 'eq.PW.png'
set style fill solid
set multiplot layout 1,1 title ""
#set xtics  0.5
#set xrange [0.3:3.2]
#set yrange [-600:400]
#set xlabel 'Time/ns'
#set xtics  0.3
#set ytics  20
set ylabel 'Press'
set xlabel 't'
#set ytics  200
p 'log_p' u (\$1*1.0E-6):(\$12*${atm2mpa}*${fG}) w lp pt 0 lw 2 lc 1 lt 1  title 'p01',\
'log_p' u (\$1*1.0E-6):(\$13*${atm2mpa}*${fG}) w lp pt 0 lw 2 lc 2 lt 1  title 'p02',\
'log_p' u (\$1*1.0E-6):(\$14*${atm2mpa}*${fG}) w lp pt 0 lw 2 lc 3 lt 1  title 'p03'


EOF
gnuplot plot.plt
rm plot.plt



cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1
#set title "Time ${t} ps"
set output 'Equalibrium.info.pp.acm.png'
set style fill solid
set multiplot layout 1,1 title ""
#set xtics  0.5
set xlabel 'Time [ns]'
set ylabel 'Press [MPa]'
#set xtics  5
#set ytics  0.005

p 'pp.10.200.dat' u (\$1*1.0E-6):(\$2*${atm2mpa}*${fG}) w lp pt 0 lw 2 lc 1 lt 1  title 'pw_x',\
'pp.10.200.acm.dat' u (\$1*1.0E-6):(\$2*${atm2mpa}*${fG}) w lp pt 0 lw 2 lc 1 lt 1  title '',\
'pp.10.200.dat' u (\$1*1.0E-6):(\$3*${atm2mpa}*${fG}) w lp pt 0 lw 2 lc 2 lt 1  title 'pw_y',\
'pp.10.200.acm.dat' u (\$1*1.0E-6):(\$3*${atm2mpa}*${fG}) w lp pt 0 lw 2 lc 2 lt 1  title '',\
'pp.10.200.dat' u (\$1*1.0E-6):(\$4*${atm2mpa}*${fG}) w lp pt 0 lw 2 lc 3 lt 1  title 'pw_z',\
'pp.10.200.acm.dat' u (\$1*1.0E-6):(\$4*${atm2mpa}*${fG}) w lp pt 0 lw 2 lc 3 lt 1  title ''

EOF
gnuplot plot.plt
rm plot.plt



#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1
#set title "Time ${t} ps"
set output 'eq.FW.png'
set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set xtics 250
#set ytics 10
#set xtics  0.5
#set xrange [0.3:3.2]
#set yrange [-600:400]
#set xlabel 'Time/ns'
#set xtics  0.3
#set ytics  20
set ylabel 'Force'
set xlabel 't'
#set ytics  200

p 'log_p' u (\$1*1.0E-6):(\$11) w lp pt 0 lw 2 lc 1 lt 1  title 'f1',\
'log_p' u (\$1*1.0E-6):(\$12) w lp pt 0 lw 2 lc 2 lt 1  title 'f2',\
'log_p' u (\$1*1.0E-6):(\$13) w lp pt 0 lw 2 lc 3 lt 1  title 'f3'

EOF
gnuplot plot.plt
rm plot.plt


