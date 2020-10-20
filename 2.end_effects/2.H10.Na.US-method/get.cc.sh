###
#####
#for AT in -0.10  0.00  0.10  0.20  0.30  0.40  0.50  0.60  0.70  0.80
for AT in -0.20 -0.30  -0.40  -0.50  -0.60  -0.70  -0.80 -0.90 -1.00
do 
cd sup.run-${AT}
tt=$(awk '{print NR}' ./COLVAR.${AT} | tail -n1 )
echo Now at sample_${AT}, length = ${tt}
awk "NR>=3 && NR<=${tt} {print}" ./COLVAR.${AT} > temp
cat COLVAR.f temp > COLVAR.10ns.${AT}
cp COLVAR.10ns.${AT} ../
rm temp
cd ..
done 
ls -l COLVAR.10ns.*
