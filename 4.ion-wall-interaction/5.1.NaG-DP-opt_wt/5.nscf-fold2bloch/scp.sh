file="/home/xuzp/WORK/ZHOUKE/5.graphene.water/1.Na/5.1.NaG-DP-opt_wt/5.nscf-fold2bloch"
scp -r xuzp@192.170.1.2:${file}/CONTCAR  ./
scp -r xuzp@192.170.1.2:${file}/OUTCAR   ./
scp -r xuzp@192.170.1.2:${file}/log      ./
scp -r xuzp@192.170.1.2:${file}/WAVECAR_spin*   ./
# scp -r xuzp@192.170.1.2:${file}/DOSCAR ./
# scp -r xuzp@192.170.1.2:${file}/AECCAR0  ./
# scp -r xuzp@192.170.1.2:${file}/AECCAR2  ./
# scp -r xuzp@192.170.1.2:${file}/CHGCAR   ./