## fit
cat << EOF > fit.plt
# f(x)=c1*exp(-x/a)+(1-c1)*exp(-x/b)
f(x)=c1*exp(-x*a*1000)+(1-c1)*exp(-x*b*100)
fit f(x) 'Cap.ave.dat' u (\$1):2 via a,b,c1

# fit f(x) 'His.Lt.10ps.dat' u (\$3):4 via a,b


EOF
gnuplot<fit.plt

rm fit.plt
