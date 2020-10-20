## fit
cat << EOF > fit.plt
f(x)=c1*exp(-x/a)+(1-c1)*exp(-x/b)
fit f(x) '01_norm_rtcf.data' u (\$1):2 via a,b,c1

#f(x)=exp(-x/a)
#fit f(x) '01_norm_rtcf.data' u (\$1):2 via a,b


EOF
gnuplot<fit.plt

rm fit.plt
