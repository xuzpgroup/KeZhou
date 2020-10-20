##
l=4
for i in 10 20 30 40 50 60 70 80 90 
do 
cd rerun.GO.bck.${i}/Str_sample.Debug.Final.GO
#cat shell.life.info.dat | tail -n +${l} | head -n ${l} 
#awk "NR==$l {print}" shell.life.info.dat
cat shell.life.info.dat
echo "----------------------"

#awk "NR>=$st && NR<=$ed {print}" log >log_press
cd ../..

done 
