/// Math function 
/// By Ke Zhou 
// ###########################################
// ###########################################
// ###########################################
#include <math.h>


int MAX_INT(int x, int y);

int MIN_INT(int x, int y);

double MAX_D(double x, double y);

double MIN_D(double x, double y);

void swapint(int *a, int *b);

void swapdouble(double *a, double *b);

double Distant(double *a, double *b);

double Distant2D(double *a, double *b);

void Vec_Minus(double *a, double *b, double *s);

void Vec_Minus_nromal(double *a, double *b, double *s);

double pbc(double ddx, double lx);

void Vec_Minus_PBC(double *a, double *b, double *s, int normlized, double *box);

double DotProduct(double *a, double *b);

int ClockwiseOrNot(double *a, double *b, double *c);

int equal_double(double x, double y, double delt);

double Distant_PBC(double *a, double *b, double *cell);

double Distant_PBC_2D(double *a, double *b, double *cell);

int BOOL(int a);

int CHAIR2ID_math(char element);

double CHAIR2MASS_math(char element);


