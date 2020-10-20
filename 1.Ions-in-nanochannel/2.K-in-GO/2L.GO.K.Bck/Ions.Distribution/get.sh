rm gz.dat
touch gz.dat
for i in 10 20 30 40 ##50 60 70 80 90 
do 
cat  gz.dat ../rerun.GO.bck.${i}/2_MSD_for_M/gas.z.dat >temp
mv temp gz.dat
done 
