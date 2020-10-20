#####################################
##### Makefile for GO generation#####
#Ke Zhou 2017.04.18 #################

main: main.o math_zhou.o radial_distribution_function.o \
	IO.o distribution.o charge.o
	gcc -o HydroDy main.o math_zhou.o \
	radial_distribution_function.o distribution.o charge.o IO.o -lm
main.o: main.c math_zhou.h radial_distribution_function.h \
	distribution.h IO.h charge.h
	gcc -c -lm main.c
math_zhou.o:  math_zhou.c  math_zhou.h
	gcc -c -lm math_zhou.c
radial_distribution_function.o: radial_distribution_function.c radial_distribution_function.h
	gcc -c -lm radial_distribution_function.c
distribution.o: distribution.c distribution.h
	gcc -c -lm distribution.c
charge.o:charge.c charge.h
	gcc -c -lm charge.c
IO.o: IO.c IO.h
	gcc -c -lm IO.c	
clean:
	rm -rf *.o
