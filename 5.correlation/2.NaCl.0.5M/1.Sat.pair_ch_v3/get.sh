# gcc -o pair pair.c  -lm  -lfftw3  -L/home/xuzp/bin/lib/fftw3/lib -I/home/xuzp/bin/lib/fftw3/include
gcc -o pair pair.c -lm -lfftw3 -L/opt/fftw3.3.3/lib -I/opt/fftw3.3.3/include
ulimit -c unlimited
./pair
rm pair
