// #####################################################################
// #                                                                   #
// # This program is used for analysing the hydration shell under nano #
// # confinement, including the statical structure and dynamical       #
// # infromation.                                                      #
// #                                                                   #
// #####################################################################

//  Program written by Zhou Ke, zhoukejustin@163.com
//  Last update 10/12/2019
//  Update: 10/12/2019, revise for entry problem 

//  Reference: ######## 



#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "math_zhou.h"
#include "radial_distribution_function.h"
#include "distribution.h"
#include "IO.h"
#include "charge.h"


#define DEBUG_1
#define MAXNATOM 10000
#define MAXWT 2000
#define Rcut 3.155
#define Rcut2nd 6.000
#define HisItems 10
#define MAXBINS 500
#define MAXT 500000
#define MAXWALLS 8000


int main()
{
	int len;
	char header_char[300];
	char filename[300];
	char filenameq[300];
	double box[3];

	box[0]=12.78;
	box[1]=12.298;
	box[2]=20.0;

	write_hander();
	sprintf(filename, "../G-WT-G-WT-1.xyz-pos-1.xyz");
	// radial_distribution_functions(120, box, filename, 8, 8);
	// radial_distribution_functions(120, box, filename, 1, 8);
	
	distributions_density(120, box, filename);

	distributions_pair_distance(120, box, filename);
	
	sprintf(filenameq, "../G-WT-G-WT--1.mulliken");
	
	distributions_charge(120, box, filenameq, filename);
	dump_into_box(box, filename);

	write_tail();
	return 0;
}



