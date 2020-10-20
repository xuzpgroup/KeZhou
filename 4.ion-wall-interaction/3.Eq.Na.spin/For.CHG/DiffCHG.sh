####
echo "A+B"
perl chgsum.pl ../For.CHG.G/CHGCAR ../For.CHG.Na/CHGCAR   
echo "AB-(A+B)"
perl chgdiff.pl ../CHGCAR ./CHGCAR_sum  
mv CHGCAR_diff CHGCAR
rm CHGCAR_sum  