#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set term png  enhanced 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Equalibrium.info.RDF.1.png'
set style fill solid
set multiplot layout 2,1 title ""
set xtics  2.0
set xlabel 'r [{\305}]'
set ylabel 'g(r)'
set xrange [0:12]
#set yrange [0:20]
set grid 

#set xtics  5
#set ytics  0.005
#set ytics 0.1

#set ytics 0.1
p './RDF/rdf.2000000.data' u (\$2):(\$3) w l lw 2 lc 1 lt 1  title 't=2.0 ns',\
'./RDF/rdf.5000000.data' u (\$2):(\$3) w l lw 2 lc 4 lt 1  title 't=5.0 ns',\
'./RDF/rdf.10000000.data' u (\$2):(\$3) w l lw 2 lc 4 lt 1  title 't=10.0 ns',\


EOF
gnuplot plot.plt
rm plot.plt


