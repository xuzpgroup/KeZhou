#!/bin/bash
## plot figure 
KJ2eV=0.0103643
KJmol2kbT=2.493
MidP=24.75
cat << EOF > plot.plt
# set term post enhanced 20 color eps
# set encoding iso_8859_1
# set term png  enhanced 
set terminal postscript eps enhanced color  font 'Arial,32' linewidth 1
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'OH.ang.eps'
set style fill solid
set multiplot layout 1,1 title ""
set size ratio 0.4
load "mycol.pal"


#set xlabel '{/Times-Italic d} [{\305}]'
set xlabel '{/Symbol q}'
set xtics  30
set ytics  0.01
# set mxtics 4
# set mytics 2
# set xrange [0.0:1.0]
# set yrange [:6]

# unset ylabel
# set lmargin 2
# set rmargin 2
# set ylabel 'F [KJ/mol]'
set ylabel '{/Arial-Italic P}({/Symbol q})'
#set ytics  0.4
#unset xlabel
# unset ytics
#set key center at 0.55, 0.8  ## set legend position
#set logscale y 
# unset label
set key font ",25"
# set key l

set style fill solid 0.5

# set label 1 at 7, 5.2 'Bulk SPC/E, our = 2.75'
# set label 2 at 7, 5.5 'Bulk TIP3P, PCCP = 6.22'
# set arrow  from 6.0, 2.75 to 24.0, 2.75 nohead lt 2 lc rgb '#DC143C'  ## ad a line

#p 'fes.dat' u (\$1*10):(\$2/${KJmol2kbT}*0.0103643) w lp pt 4 ps 2 lw 2 lc 1 lt 1  title 'd',\

p 'O-H.ang.dat' u (\$1):(\$2*-1) \
w l ls 11 lw 8 title 'Bottom',\
'O-H.ang.dat' u (\$1):(\$3*1) \
w l ls 12 lw 8   title 'Top',\


### smooth  sbezier

EOF
gnuplot plot.plt
rm plot.plt

cat << EOF > plot.plt
# set term png  enhanced 
set terminal postscript eps enhanced color  font 'Arial,25' linewidth 1
set encoding iso_8859_1
set output 'OH.ang.2D.eps'

set style fill solid
set multiplot layout 1,1 title "{/Arial-Italic E} = 0 V/nm"
set size ratio 0.4

#set format x "%.e "
set xtics 30
# set ytics 0.2
# set xtics 0.2
set xrange [0:180]
set yrange [0:10]
# set xlabel '{/Symbol q}'
set xlabel 'ang'
set ylabel '{/Arial-Italic z} ({\305})'

#set ytics 0.1
# load 'Moreland.plt'
# load 'plasma.pal'
# load 'parula.pal'
# load 'viridis.pal'
load 'jet.palette'


# unset ylabel
# unset ytics
# unset cbtics
set cbrange [0:0.003]
set cblabel "Count"

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


# set style fill transparent solid 0.5 ##noborder

# w filledcurves fs transparent solid 0.6 \


p 'O-H.ang.z.dat' u (\$1/1):(\$2/1):(\$3)  w image ti '',\


EOF
gnuplot plot.plt
rm plot.plt



# for i in *.eps 
# do \
echo "Convert to png"
# convert -density 500 Density.profile.AIMD.eps Density.profile.AIMD.png 
# done 





