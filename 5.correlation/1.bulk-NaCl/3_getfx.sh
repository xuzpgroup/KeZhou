#!/bin/bash
#awk  '/Step Temp TotEng Lx Ly Lz Volume Press/{print NR}' log >line.1
sed -n '/Step /=' log > line.1
sed -n '/Loop time/=' log > line.2
a1=$(sed -n '8p' line.1)
b1=$(sed -n '8p' line.2)
b1=$(awk '{print NR}' log | tail -n1 )
echo ${a1} 
echo ${b1}  
st=`expr $a1 + 0`
ed=`expr $b1 - 1`
awk "NR>=$st && NR<=$ed {print}" log >log_p
rm line.1 line.2
#a1=$(sed -n '8p' line.1)
#b1=$(sed -n '8p' line.2)
#echo ${a1} 
#echo ${b1}  
#st=`expr $a1 + 2`
#ed=`expr $b1 - 1`
#awk "NR>=$st && NR<=$ed {print}" log >log_2

#cat log_1 log_2 >log_press_bak
#rm log_1 log_2

#tt=$(awk '{print NR}' log_press | tail -n1 )

#nn=10000 ## statics nn last lines
#rm -rf log_press_only
#rm -rf log_cell_only
#echo ${tt} >log_press_only
#echo ${nn} >>log_press_only
#echo ${tt} >log_cell_only
#@echo ${nn} >>log_cell_only

#awk '{printf $1 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\n" }' log_press_bak >log_msd_1_only
#awk '{printf $1 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\n" }' log_press_bak >log_msd_2_only
#awk '{printf $2 "\t" $3 "\n"}' log_press >>log_cell_only

# calculate press 
#gcc -o 1_cal_press.x 1_cal_press.cpp
#./1_cal_press.x
# plot press with time
#sh 1_plot_press.sh




#echo "#############"

