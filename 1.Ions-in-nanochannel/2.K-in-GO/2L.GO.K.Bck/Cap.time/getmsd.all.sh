## CH4
p1="../rerun.GO.bck.10/Str_sample.Debug.Final.GO/Capture_time/01_norm_rtcf.data"
p2="../rerun.GO.bck.20/Str_sample.Debug.Final.GO/Capture_time/01_norm_rtcf.data"
p3="../rerun.GO.bck.30/Str_sample.Debug.Final.GO/Capture_time/01_norm_rtcf.data"
p4="../rerun.GO.bck.40/Str_sample.Debug.Final.GO/Capture_time/01_norm_rtcf.data"
p5="../rerun.GO.bck.50/Str_sample.Debug.Final.GO/Capture_time/01_norm_rtcf.data"
p6="../rerun.GO.bck.60/Str_sample.Debug.Final.GO/Capture_time/01_norm_rtcf.data"
p7="../rerun.GO.bck.70/Str_sample.Debug.Final.GO/Capture_time/01_norm_rtcf.data"
p8="../rerun.GO.bck.80/Str_sample.Debug.Final.GO/Capture_time/01_norm_rtcf.data"
p9="../rerun.GO.bck.90/Str_sample.Debug.Final.GO/Capture_time/01_norm_rtcf.data"

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
 msd.case.6.dat msd.case.7.dat msd.case.8.dat msd.case.9.dat > 01_norm_rtcf.all.data
rm msd.case.1.dat msd.case.2.dat msd.case.3.dat msd.case.4.dat msd.case.5.dat 
rm msd.case.6.dat msd.case.7.dat msd.case.8.dat msd.case.9.dat 
