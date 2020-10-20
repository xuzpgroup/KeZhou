cat << EOF > plot.plt
set term png  enhanced  font "arial,15" fontscale 1.0  ##size 512, 280 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'His.Lt.all.dx.png'
set style fill solid
set multiplot layout 1,1 title "Probility distribution of particle dispalcement"

# set lmargin 8
# set rmargin 8
# set tmargin 8
# set bmargin 3


set ylabel 'Distribution'
set xlabel 'Diffusion length [{\305}]'
#unset ytics
set xrange [0.0:10.0]
#set yrange [0:1.2]
#set ytics  0.2
#set logscale y
#set yrange [0.01:2]
p 'His.Lt.10ps.dat' u (\$1):(\$2) w p pt 7 ps 1 lw 1 lc 1 lt 1  title 'dx, dt=10 ps',\
1.2766*exp(-x*x/0.552042) w l lt 1 lc 1 lw 2  title '',\
'./20ps/His.Lt.20ps.dat' u (\$1):(\$2) w p pt 7 ps 1 lw 1 lc 2 lt 1  title 'dt=20 ps',\
0.992096*exp(-x*x/0.810267) w l lt 1 lc 2 lw 2  title '',\
'./30ps/His.Lt.30ps.dat' u (\$1):(\$2) w p pt 7 ps 1 lw 1 lc 3 lt 1  title 'dt=30 ps',\
0.831686*exp(-x*x/1.09077) w l lt 1 lc 3 lw 2  title '',\
'./40ps/His.Lt.40ps.dat' u (\$1):(\$2) w p pt 7 ps 1 lw 1 lc 4 lt 1  title 'dt=40 ps',\
0.682674*exp(-x*x/1.72278) w l lt 1 lc 4 lw 2  title '',\
'./50ps/His.Lt.50ps.dat' u (\$1):(\$2) w p pt 7 ps 1 lw 1 lc 5 lt 1  title 'dt=50 ps',\
0.581574*exp(-x*x/2.52874) w l lt 1 lc 5 lw 2  title '',\
'./100ps/His.Lt.100ps.dat' u (\$1):(\$2) w p pt 7 ps 1 lw 1 lc 6 lt 1  title 'dt=100 ps',\
0.386349*exp(-x*x/6.05985) w l lt 1 lc 6 lw 2  title ''


# 'His.Lt.10ps.dat' u (\$3):(\$4) w p pt 7 ps 1 lw 1 lc 2 lt 1  title 'dy',\
# 0.516531*exp(-x*x/4.74784) w l lt 1 lc 2 lw 2  title 'Gauss fit'



EOF
gnuplot plot.plt
rm plot.plt


cat << EOF > plot.plt
set term png  enhanced  font "arial,15" fontscale 1.0  ##size 512, 280 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'His.Lt.all.dy.png'
set style fill solid
set multiplot layout 1,1 title "Probility distribution of particle dispalcement"

# set lmargin 8
# set rmargin 8
# set tmargin 8
# set bmargin 3


set ylabel 'Distribution'
set xlabel 'Diffusion length [{\305}]'
#unset ytics
set xrange [0.0:10.0]
#set yrange [0:1.2]
set ytics  0.2
p 'His.Lt.10ps.dat' u (\$3):(\$4) w p pt 7 ps 1 lw 1 lc 1 lt 1  title 'dy, dt=10 ps',\
0.944736*exp(-x*x/1.27532) w l lt 1 lc 1 lw 2  title '',\
'./20ps/His.Lt.20ps.dat' u (\$3):(\$4) w p pt 7 ps 1 lw 1 lc 2 lt 1  title 'dt=20 ps',\
0.710915*exp(-x*x/2.23421) w l lt 1 lc 2 lw 2  title '',\
'./30ps/His.Lt.30ps.dat' u (\$3):(\$4) w p pt 7 ps 1 lw 1 lc 3 lt 1  title 'dt=30 ps',\
0.605383*exp(-x*x/3.02184) w l lt 1 lc 3 lw 2  title '',\
'./40ps/His.Lt.40ps.dat' u (\$3):(\$4) w p pt 7 ps 1 lw 1 lc 4 lt 1  title 'dt=40 ps',\
0.546236*exp(-x*x/3.56234) w l lt 1 lc 4 lw 2  title '',\
'./50ps/His.Lt.50ps.dat' u (\$3):(\$4) w p pt 7 ps 1 lw 1 lc 5 lt 1  title 'dt=50 ps',\
0.506916*exp(-x*x/4.01497) w l lt 1 lc 5 lw 2  title '',\
'./100ps/His.Lt.100ps.dat' u (\$3):(\$4) w p pt 7 ps 1 lw 1 lc 6 lt 1  title 'dt=100 ps',\
0.371065*exp(-x*x/7.00467) w l lt 1 lc 6 lw 2  title ''


# 'His.Lt.10ps.dat' u (\$3):(\$4) w p pt 7 ps 1 lw 1 lc 2 lt 1  title 'dy',\
# 0.516531*exp(-x*x/4.74784) w l lt 1 lc 2 lw 2  title 'Gauss fit'



EOF
gnuplot plot.plt
rm plot.plt

