gcc -o pair pair.c  -lm  -lfftw3 -L/opt/fftw3.3.3/lib -I/opt/fftw3.3.3/include
./pair
rm pair
