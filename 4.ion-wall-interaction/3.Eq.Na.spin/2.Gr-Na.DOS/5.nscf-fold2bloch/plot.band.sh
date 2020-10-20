#!/bin/bash
## plot figure 
cat << EOF > plot.plt
set terminal postscript eps enhanced color  font 'Arial,20' linewidth 1 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'band.sp.eps'
set style fill solid
set multiplot layout 1,2 title ""
set size ratio 1.8

# load "mycol.pal"
# load "Moreland.plt"
# load "viridis.pal"
set lmargin 8 
set rmargin 2 
# set key font ',20'
# set key b
# set cbrange [0:4]


# set cbrange [0.0:2.0]
set yrange [-10:5]
set ytics  5
# set xtics  1
set ylabel '{/Arial-Italic E} (eV)'
unset xlabel
unset xtics
set xtics ("G" 0.0, "M"  0.2339, "K" 0.3689, "G" 0.6390)
  
    

# set pm3d map
set pm3d interpolate 0,0

# set palette defined (0 "black", 1 "yellow", 2 "red", 3 "blue")
# set palette defined (0 "green", 1 "blue", 2 "red", 3 "orange")
# set palette defined (0 "blue", 1 "green", 2 "yellow", 3 "red")
# set palette defined (0 "#16469D", 1 "#3bb9d7", 2 "yellow", 3 "red")
load "jet.palette"


plot "WAVECAR_spin1.f2b.bin" binary matrix with image title "",\
     0 title "" lc "white" lw 2 lt 2 dt 2,\
     "HKP.dat" u 1:2 w l lw 4 dt 3 lc "white" ti "" 


plot "WAVECAR_spin2.f2b.bin" binary matrix with image title "",\
     0 title "" lc "white" lw 2 lt 2 dt 2,\
     "HKP.dat" u 1:2 w l lw 4 dt 3 lc "white" ti "" 

EOF
gnuplot plot.plt
rm plot.plt

file="band.sp"
convert -density 500 ${file}.eps ${file}.Na@G-E0-.png 