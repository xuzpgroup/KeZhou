## CH4
p1="../rerun.GO.bck.10/Velocity_accoralation/03.trap.acf.dat"
p2="../rerun.GO.bck.20/Velocity_accoralation/03.trap.acf.dat"
p3="../rerun.GO.bck.30/Velocity_accoralation/03.trap.acf.dat"
p4="../rerun.GO.bck.40/Velocity_accoralation/03.trap.acf.dat"
p5="../rerun.GO.bck.50/Velocity_accoralation/03.trap.acf.dat"
p6="../rerun.GO.bck.60/Velocity_accoralation/03.trap.acf.dat"
p7="../rerun.GO.bck.70/Velocity_accoralation/03.trap.acf.dat"
p8="../rerun.GO.bck.80/Velocity_accoralation/03.trap.acf.dat"
p9="../rerun.GO.bck.90/Velocity_accoralation/03.trap.acf.dat"

awk '{printf $1 " " $2 " " "\n"  }'  ${p1} >msd.case.1.dat
awk '{printf $2 " " "\n"  }'  ${p2} >msd.case.2.dat
awk '{printf $2 " " "\n"  }'  ${p3} >msd.case.3.dat
awk '{printf $2 " " "\n"  }'  ${p4} >msd.case.4.dat
awk '{printf $2 " " "\n"  }'  ${p5} >msd.case.5.dat
awk '{printf $2 " " "\n"  }'  ${p6} >msd.case.6.dat
awk '{printf $2 " " "\n"  }'  ${p7} >msd.case.7.dat
awk '{printf $2 " " "\n"  }'  ${p8} >msd.case.8.dat
awk '{printf $2 " " "\n"  }'  ${p9} >msd.case.9.dat

paste msd.case.1.dat msd.case.2.dat msd.case.3.dat msd.case.4.dat msd.case.5.dat \
 msd.case.6.dat msd.case.7.dat msd.case.8.dat msd.case.9.dat > msd.all.data
rm msd.case.1.dat msd.case.2.dat msd.case.3.dat msd.case.4.dat msd.case.5.dat 
rm msd.case.6.dat msd.case.7.dat msd.case.8.dat msd.case.9.dat 
