## fit
cat << EOF > fit.plt
f(x)=a*exp(-x*x/b)
fit f(x) 'His.Lt.20ps.dat' u (\$1):2 via a,b
fit f(x) 'His.Lt.20ps.dat' u (\$3):4 via a,b


EOF
gnuplot<fit.plt

rm fit.plt
