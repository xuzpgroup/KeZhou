tac OUTCAR | grep -m 1 -B 0 "free  energy   TOTEN" | awk '{print $5}'
