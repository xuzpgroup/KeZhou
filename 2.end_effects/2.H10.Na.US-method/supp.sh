#####
for AT in -0.10  0.00  0.10  0.20  0.30  0.40  0.50  0.60  0.70  0.80
do 
cd sample_${AT}
mkdir sup.run 
cp t.v.restart ./sup.run/re.restart
cp ../re.inp ./sup.run/eq.GO.ions.in
cp FF.water.ion.prm tsjobs plumed.dat ./sup.run
cd sup.run
bsub<tsjobs 
cd ../
cd ..
done 


