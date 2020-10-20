###
for i in 10 20 30 40 50 60 70 80 90 
do 
mkdir rerun.GO.bck.${i}
sed -e 's/YY/'${i}'/g' ./K.temp.1 > ./rerun.GO.bck.${i}/G2.K.data
sed -e 's/XX/'${i}'/g' ./tsing-job > ./rerun.GO.bck.${i}/tsing-job
cp tes  eq.GO.ions.in ./rerun.GO.bck.${i}
cd rerun.GO.bck.${i}

cd ..

done
