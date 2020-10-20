#!/bin/bash
## plot figure 
H=10.0
SHFT=25.5443
ST=52.4567
fes="../free.out"
cat << EOF > plot.plt
# set term png  enhanced 
set terminal postscript eps enhanced color  font 'Arial,20' linewidth 1
set encoding iso_8859_1

#set title "Time ${t} ps"
set output '2D.his.ions.eps'

set style fill solid
set multiplot layout 1,1 title ""

#set format x "%.e "
#set xtics 250
set ytics 0.4
set xtics 0.5
set xrange [-1:1]
set yrange [0:1.36]
set ylabel '{/Arial-Italic z} (nm)'
set xlabel '{/Arial-Italic x} (nm)'

#set ytics 0.1
# load 'Moreland.plt'
# load 'plasma.pal'
load 'parula.pal'
# load 'viridis.pal'

# unset ylabel
# unset ytics
# set cbrange [0:5.0]
set cblabel "Count"

# set xrange [0.0-${ST}/10:8.400-${ST}/10]
# set yrange [0.0:2.00]

set pm3d map
set view map scale 1
# set samples 25, 25
# set isosamples 26, 26
			### fro coutour lines
# unset surface 
# set contour base
# # set cntrlabel  format '%8.3g' font ',7' start 2 interval 20
# # set cntrparam order 8
# # set cntrparam bspline
# # set cntrparam levels auto 10
# set style increment user
# do for [i=1:7] {set style line i lc rgb "black"}
# set cntrparam level incremental 0.0, 1.0, 3.0


set size ratio 0.3333333

# set key out 
# set key right
# set key at 0,0 
# set key box opaque
# set key width 2
# set key height 0.5

# set x2range [0.0:84.00]
# set y2range [0:30.0]
# set y2label 'PMF [KCal/mol]'
# set y2tics 5


set style fill transparent solid 0.5 noborder

# w filledcurves fs transparent solid 0.6 \


p 'His.2D.ions.dat' u ((\$1)/10-5.24924):(\$2/10):(\$3)  w image ti '',\
'./His.1D.dat' u ((\$1)/10-5.24924):(\$2/5000) \
w l lw 4 lc rgb 'red' lt 1  title '',\
'../1.H6.8.M1.Na.temp.pmf.2nd.data' u ((\$1)/10-1.0):(\$2/4+0.2)\
  w  l lw 3 lc 3 ti '',\

# 'His.2D.dat' u ((\$4-${ST}+84.0)/10):(\$5/10):(\$6/1000)  w image ti '',\
# 'oxided.xyz' u (\$1):(\$2)  w p pt 7 ps 2 lc rgb 'red' ti ''

# 'Water.distribution.data' u (\$3*50+10.0):(\$1+${H}*0.0) w l lw 2 lc rgb "white" lt 1 ti '',\
# 'Water.distribution.data' u (\$3*50+10.0):(\$1+${H}*2.0) w l lw 2 lc rgb "white" lt 1 ti '',\
# 'Water.distribution.data' u (\$3*50+10.0):(\$1+${H}*4.0) w l lw 2 lc rgb "white" lt 1 ti '',\

# '../free.out' u (\$1*10+32.5):(\$2) axes x2y2 w l lw 3 lc 3 ti 'Na^+',\

# splot 'His.2D.dat' u (\$4):(\$5):(\$6/1000)  with lines lc 0 lw 3 ti '',\


# set pm3d map
# sp 'His.2D.dat' u (\$4):(\$5):(\$6/1000),\

# # 'wall.xyz' u (\$2):(\$4)  w p pt 6 ps 4 lc 2 ti '',\


EOF
gnuplot plot.plt
rm plot.plt



# cat << EOF > plot.plt
# # set term png  enhanced 
# set terminal postscript eps enhanced color  font 'Arial,20' linewidth 1
# set encoding iso_8859_1

# #set title "Time ${t} ps"
# set output '1D.his.z.eps'

# set style fill solid
# set multiplot layout 1,3 title ""

# set xlabel 'Density [kg/m^3]'
# set ylabel 'Z [{\305}]'
# set xtics  1
# set ytics  5

# # load 'Moreland.plt'
# load 'plasma.pal'
# # load 'parula.pal'
# # load 'viridis.pal'

# # set xrange [3.0:0.0]
# # set yrange [-15.0:15.000]


# # p 'His.2D.dat' u (\$1):(\$2):(\$3) w image
# p 'His.1D.dat' u (\$4):(\$3-15.0) w l lc 0 lw 2 title ''



# EOF
# gnuplot plot.plt
# rm plot.plt



