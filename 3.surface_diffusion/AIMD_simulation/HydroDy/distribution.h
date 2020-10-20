/// Functional tools 
/// By Ke Zhou 
// ###########################################
// ###########################################
// ###########################################

#include "math_zhou.h"
#include <math.h>

#define MAX_ATOMS 300
#define PI 3.1415926


double CHAIR2MASS_disp(char element);

int CHAIR2ID_disp(char element);

void HistogramIt_1D(double (*his)[3], double samples, double st, \
	double dt, double fac, int ind);

void distributions_density(int nwalls, \
    double *box, char *filename);

void distributions_pair_distance(int nwalls, \
    double *box, char *filename);

