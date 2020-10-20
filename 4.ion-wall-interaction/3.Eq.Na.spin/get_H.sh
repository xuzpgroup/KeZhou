# gcc -o chg_atom chg_atom.c 
# ./chg_atom 
# rm chg_atom


echo "###############"
gcc -o cal_height cal_height.c -lm
./cal_height
rm cal_height

# echo "###############"
# tail -n 8 OSZICAR

echo "####### Mu ######"
grep "mag=" log  | tail -n1 | awk '{printf $10 "\n" }'

echo "####### Dip ######"
grep "dipolmoment" OUTCAR | tail -n1 | awk '{printf $4 "\n" }' 


