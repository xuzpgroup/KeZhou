#!/bin/bash
####
gas="rerun"
#rm Bind.total.${gas}.dat
#echo Bind.total.${gas}.dat
workpath="/home/zhouke/WORK_FILE/17-04-11-ION-Transport/4.2_CMD_run/1.1.GO.6.8.Na/00.3.Browian_anas"
cd ../
f1="Lt.20.dat"
f2="Lt.40.dat"


cp 2_MSD_for_gas/${f1} ./00.3.Browian_anas/Lt.tot.${f1}
cp 2_MSD_for_gas/${f2} ./00.3.Browian_anas/Lt.tot.${f2}


for ((i=5;i<=45;i=$i+5))
do 
	echo "Now ${gas}.${i}" 
	cp -r ${workpath}/3_Brownian.sh ./${gas}.${i}/2_MSD_for_gas
	cp -r ${workpath}/Borwnian.c    ./${gas}.${i}/2_MSD_for_gas
	cd ${gas}.${i}/2_MSD_for_gas
	sh 3_Brownian.sh
	cat ./${f1} ../../00.3.Browian_anas/Lt.tot.${f1} > ../../00.3.Browian_anas/Temp.1
	mv ../../00.3.Browian_anas/Temp.1 ../../00.3.Browian_anas/Lt.tot.${f1}
	cat ./${f2} ../../00.3.Browian_anas/Lt.tot.${f2} > ../../00.3.Browian_anas/Temp.1
	mv ../../00.3.Browian_anas/Temp.1 ../../00.3.Browian_anas/Lt.tot.${f2}
	cd ../../
done 
