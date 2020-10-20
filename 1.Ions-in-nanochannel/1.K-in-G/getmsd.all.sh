## CH4
p1="2_MSD_for_gas/msd.data"
p2="rerun.5/2_MSD_for_gas/msd.data"
p3="rerun.10/2_MSD_for_gas/msd.data"
p4="rerun.15/2_MSD_for_gas/msd.data"
p5="rerun.20/2_MSD_for_gas/msd.data"
p6="rerun.25/2_MSD_for_gas/msd.data"
p7="rerun.30/2_MSD_for_gas/msd.data"
p8="rerun.35/2_MSD_for_gas/msd.data"
p9="rerun.40/2_MSD_for_gas/msd.data"
p10="rerun.45/2_MSD_for_gas/msd.data"

awk '{printf $1 " " $2 " " "\n"  }'  ${p1} >msd.case.1.dat
awk '{printf $2 " " "\n"  }'  ${p2} >msd.case.2.dat
awk '{printf $2 " " "\n"  }'  ${p3} >msd.case.3.dat
awk '{printf $2 " " "\n"  }'  ${p4} >msd.case.4.dat
awk '{printf $2 " " "\n"  }'  ${p5} >msd.case.5.dat
awk '{printf $2 " " "\n"  }'  ${p6} >msd.case.6.dat
awk '{printf $2 " " "\n"  }'  ${p7} >msd.case.7.dat
awk '{printf $2 " " "\n"  }'  ${p8} >msd.case.8.dat
awk '{printf $2 " " "\n"  }'  ${p9} >msd.case.9.dat
awk '{printf $2 " " "\n"  }'  ${p10} >msd.case.10.dat

paste msd.case.1.dat msd.case.2.dat msd.case.3.dat msd.case.4.dat msd.case.5.dat \
 msd.case.6.dat msd.case.7.dat msd.case.8.dat msd.case.9.dat msd.case.10.dat > msd.all.data
rm msd.case.1.dat msd.case.2.dat msd.case.3.dat msd.case.4.dat msd.case.5.dat 
rm msd.case.6.dat msd.case.7.dat msd.case.8.dat msd.case.9.dat msd.case.10.dat 
