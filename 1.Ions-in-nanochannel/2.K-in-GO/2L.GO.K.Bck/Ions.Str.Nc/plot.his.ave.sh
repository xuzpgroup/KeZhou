rm ave.cal
gcc -o ave.cal ave.cal.c  -lm 
./ave.cal

### plot 
#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Nc.his.ave.png'
set style fill solid
set multiplot layout 3,1 title ""


set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
#set xtic rotate by -45 scale 0

#set format x "%.e "
#set xtics 250
#set ytics 10
#set xtics  2
#set xrange [1:8]
#set yrange [3:-12]

set ylabel 'P(N_c) [%]'
#set xlabel 'Snopshot'
#set xtics  5
#set ytics  0.005
#set ytics 0.1

#set ytics 0.1
p 'His.nc.ave.dat' u (\$2*100):xtic(1)  title ''  #ti col 


### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt
