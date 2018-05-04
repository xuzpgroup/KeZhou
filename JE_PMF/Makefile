#####################################
##### Makefile for GO generation#####
#Ke Zhou 2017.04.18 #################

JE_PMF: JE_PMF.o math_zhou.o 
	gcc -o JE_PMF JE_PMF.o math_zhou.o -lm
JE_PMF.o: JE_PMF.c math_zhou.h
	gcc -c -lm JE_PMF.c
math_zhou.o:  math_zhou.c  math_zhou.h
	gcc -c -lm math_zhou.c
clean:
	rm -rf *.o
