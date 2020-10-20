#####
for AT in -0.10  ###0.00  0.10  0.20  0.30  0.40  0.50  0.60  0.70  0.80
do 
cd sample_${AT}

tt=$(awk '{print NR}' ./sup.run/COLVAR.${AT} | taATl -n1 )
awk "NR>=3 && NR<=${tt} {print}" ./sup.run/COLVAR.${AT} > temp
cat COLVAR.${AT} temp > COLVAR.10ns.${AT}
rm temp
cd ..
done 
