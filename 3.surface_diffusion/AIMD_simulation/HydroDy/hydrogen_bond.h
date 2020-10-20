/// Functional tools 
/// By Ke Zhou 
// ###########################################
// ###########################################
// ###########################################

#include "math_zhou.h"
#include <math.h>

#define MAX_ATOMS 300
#define PI 3.1415926

void HistogramIt_ang(double (*his)[2], double samples, double st, \
    double dt, double fac, int ind);

void HistogramIt_ang_2D(double (*his)[800], \
	double samples1, double st1, double dt1, \
	double samples2, double st2, double dt2, double fac);

void HB_Rotation(int nwalls, \
    double *box, char *filename);


