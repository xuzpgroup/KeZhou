ifort -o fft vacf_fft.f90    -lfftw3   
##  -lfftw3 -L/home/xuzp/bin/lib/fftw3/lib -I/home/xuzp/bin/lib/fftw3/include
echo "Don FFT ..."
./fft 
gfortran -o trap trap.f90
./trap 
#rm trap 
