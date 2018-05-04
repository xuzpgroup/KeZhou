/// Math function 
/// By Ke Zhou 
// ###########################################
// ###########################################
// ###########################################

#include "math_zhou.h"
#include <math.h>


int MAX_INT(int x, int y)
{
    if (x>y) return x;
    return y;
    // return x>y?x:b;
}
int MIN_INT(int x, int y)
{
    if (x<y) return x;
    return y;
}
double MAX_D(double x, double y)
{
    if (x>y) return x;
    return y;
}
double MIN_D(double x, double y)
{
    if (x<y) return x;
    return y;
}
void swapint(int *a, int *b)
{
    int ss;
    ss=*a;     *a=*b;    *b=ss;
}
void swapdouble(double *a, double *b)
{
    double ss;
    ss=*a;     *a=*b;    *b=ss;
}
double DABS(double x)
{
	if (x<0) return -x;
	return x;
}

double Distant(double *a, double *b)
{
    double r;
    double d1, d2, d3;
    d1=(a[0]-b[0]);
    d2=(a[1]-b[1]);
    d3=(a[2]-b[2]);
    r=sqrt(d1*d1+d2*d2+d3*d3);  
    return r;
}

double Distant2D(double *a, double *b)
{
    double r;
    double d1, d2;
    d1=(a[0]-b[0]);
    d2=(a[1]-b[1]);
    r=sqrt(d1*d1+d2*d2);  
    return r;
}

double Distant_PBC(double *a, double *b, double *cell)
{
    double r;
    double d1, d2, d3;
    d1=DABS(a[0]-b[0]); if (d1>cell[0]/2) d1=cell[0]-d1;
    d2=DABS(a[1]-b[1]); if (d2>cell[1]/2) d2=cell[1]-d2;
    d3=DABS(a[2]-b[2]); if (d3>cell[2]/2) d3=cell[2]-d3;
    r=sqrt(d1*d1+d2*d2+d3*d3);  
    return r;	
}

double Distant_PBC_exp_x(double *a, double *b, double *cell)
{
    double r;
    double d1, d2, d3;
    d1=DABS(a[0]-b[0]); // if (d1>cell[0]/2) d1=cell[0]-d1;
    d2=DABS(a[1]-b[1]); if (d2>cell[1]/2) d2=cell[1]-d2;
    d3=DABS(a[2]-b[2]); if (d3>cell[2]/2) d3=cell[2]-d3;
    r=sqrt(d1*d1+d2*d2+d3*d3);  
    return r;   
}

double Distant_PBC_2D(double *a, double *b, double *cell)
{
    double r;
    double d1, d2;
    d1=DABS(a[0]-b[0]); if (d1>cell[0]/2) d1=cell[0]-d1;
    d2=DABS(a[1]-b[1]); if (d2>cell[1]/2) d2=cell[1]-d2;
    r=sqrt(d1*d1+d2*d2);  
    return r;   
}

void Vec_Minus(double *a, double *b, double *s)
{
	s[0]=a[0]-b[0];
	s[1]=a[1]-b[1];
	s[2]=a[2]-b[2];	
}

double DotProduct(double *a, double *b)
{
	double s;
	s=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	return s;
}

int ClockwiseOrNot(double *a, double *b, double *c)
{
	double vx, vy, ux, uy;
	double k;
	vx=(b[0]-a[0]);
	vy=(b[1]-a[1]);
	ux=(c[0]-b[0]);
	uy=(c[1]-b[1]);
	k=vx*uy-ux*vy;
	if (k>0) return -1;   // counter clockwise
	return 1;
}

int equal_double(double x, double y, double delt)
{
	if (DABS(x-y)<delt) return 1;
	return 0;
}

int BOOL(int a)
{
    if (a==0) return 0;
    return 1;
}

double array_ave(int len, double *a)
{
    int i;
    double s=0;
    for (i=0;i<len;i++)
        s=s+a[i];
    s=s/len;
    return s;
}

double array_square_ave(int len, double *a)
{
    int i;
    double s=0;
    for (i=0;i<len;i++)
        s=s+a[i]*a[i];
    s=s/len;
    return s;
}

double array_bolzman_ave(int len, double *a, double kb, double temp)
{
    int i;
    double s=0;
    double beta;
    beta=1.0/kb/temp;

    for (i=0;i<len;i++)
        s=s+exp(-a[i]*beta);
    s=s/len;
    return s;
}


