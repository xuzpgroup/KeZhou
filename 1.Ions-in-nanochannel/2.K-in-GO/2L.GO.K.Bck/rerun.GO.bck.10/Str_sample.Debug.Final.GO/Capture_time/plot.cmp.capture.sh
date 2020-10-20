T=0.75
PI=3.1415926

fNa="/home/zhou/WORK_FILE/17-04-11-ION-Transport/4.4._CMD.GO/1.0_Na/2L.GO.Na.Bck/rerun.GO.bck.10/Str_sample.Debug.Final.GO/Capture_time"
fK="/home/zhou/WORK_FILE/17-04-11-ION-Transport/4.4._CMD.GO/2.0_K/2L.GO.K.Bck/rerun.GO.bck.10/Str_sample.Debug.Final.GO/Capture_time"
#!/bin/bash
## plot figure 
cat << EOF > plot.plt
#set terminal postscript eps color solid linewidth 2 "Helvetica" 26
# set terminal pngcairo  transparent enhanced font "arial,12" fontscale 1.0  ##size 512, 280 
set terminal pngcairo  enhanced font "arial,15" fontscale 1.0  ##size 512, 280 
# set term png  enhanced  font "arial,15" fontscale 1.0  ##size 512, 280 
set encoding iso_8859_1

#set title "Time ${t} ps"
set output 'Capture.RTCF.all.png'
set style fill solid
# set multiplot layout 1,1 title "Capture time correlation function"
set multiplot layout 1,1 title "Capture time: Residence time correlation function"


#set xrange [0:100]
set yrange [0:1.5]
set ytics  0.5
set xtics  25
#set yrange [0:5]
set xlabel 't [ps]'
#set ylabel 'Q ({/Symbol-Oblique w})'
set ylabel '<v(0)v(t)> '
#set xtics  5
#unset xlabel
#unset xtics
#set key center at 130, 4.0  ## set legend position
#set logscale

#p 'Q.omega.data' u (\$1*2*${PI}/${T}*1.0E-3):(\$3) w p pt 0 lw 2 lc 0 lt 1  title ''
p '${fNa}/01_norm_rtcf.data' u (\$1):(\$2) w l lw 2 lc 1 lt 1  title 'Na^+',\
 0.381748*exp(-x/716.762)+(1-0.381748)*exp(-x/39.0321) w l lt 1 dt 3 lc 1 lw 2  title "Fitting",\
'${fK}/01_norm_rtcf.data' u (\$1):(\$2) w l lw 2 lc 2 lt 1  title 'K^+',\
 0.869243*exp(-x/228.182)+(1-0.869243)*exp(-x/8.31419) w l lt 1 dt 3 lc 2 lw 2  title "Fitting"


EOF
gnuplot plot.plt
rm plot.plt


