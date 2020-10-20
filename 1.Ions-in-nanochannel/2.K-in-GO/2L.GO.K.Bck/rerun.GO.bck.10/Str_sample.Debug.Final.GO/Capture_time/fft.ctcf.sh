gcc -o ctcf fft.ctcf.c   -lfftw3 -lm  -L/opt/fftw3.3.3/lib -I/opt/fftw3.3.3/include

##Could you try -lfftw3 -lm instead of -lm -lfftw3
##
##  -lfftw3 -lm -L/home/xuzp/bin/lib/fftw3/lib -I/home/xuzp/bin/lib/fftw3/include
./ctcf
rm ctcf

