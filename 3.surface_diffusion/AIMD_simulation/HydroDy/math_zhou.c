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
    double d[3];
    int i;

    for (i=0;i<3;i++)
    {
        d[i]=a[i]-b[i];
        d[i]=cell[i]*ceil(d[i]/cell[i])-d[i];
        if (d[i]>cell[i]/2) d[i]=cell[i]-d[i];
    }    
    r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);  
    return r;	
}

double Distant_PBC_2D(double *a, double *b, double *cell)
{
    double r;
    double d[3];
    int i;

    for (i=0;i<2;i++)
    {
        d[i]=a[i]-b[i];
        d[i]=cell[i]*ceil(d[i]/cell[i])-d[i];
        if (d[i]>cell[i]/2) d[i]=cell[i]-d[i];
    }    
    r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);  
    return r;     
}

void Vec_Minus(double *a, double *b, double *s)
{
	s[0]=a[0]-b[0];
	s[1]=a[1]-b[1];
	s[2]=a[2]-b[2];	
}

void Vec_Minus_nromal(double *a, double *b, double *s)
{
    double l;
    s[0]=a[0]-b[0];
    s[1]=a[1]-b[1];
    s[2]=a[2]-b[2]; 
    l=sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
    s[0]/=l;
    s[1]/=l;
    s[2]/=l;       
}

double pbc(double ddx, double lx)
{
    double dx;
    dx=ddx;
    while (dx>0.5*lx) dx=dx-lx;
    while (dx<-0.5*lx) dx=dx+lx;
    return dx;
}


void Vec_Minus_PBC(double *a, double *b, double *s, int normlized, double *box)
{
    double l;
    // s=a-b
    s[0]=pbc(a[0]-b[0], box[0]);
    s[1]=pbc(a[1]-b[1], box[1]);
    s[2]=pbc(a[2]-b[2], box[2]);    
    l=sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
    if (normlized)
    {
        s[0]/=l;
        s[1]/=l;
        s[2]/=l;
    }
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



// ##################################
double CHAIR2MASS_math(char element)
{
    if (element=='C') return 12.00;
    if (element=='H') return 1.00;
    if (element=='O') return 16.00;   
    return -1;     
}

int CHAIR2ID_math(char element)
{
    if (element=='C') return 12;
    if (element=='H') return 1;
    if (element=='O') return 8;   
    return -1;     
}