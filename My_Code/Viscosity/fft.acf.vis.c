// calculate auto-correlation Function by FFT
// For viscosity calculation 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fftw3.h"


#define  MAX_BIN 3
#define  MAX_DIM 3
#define  MAX_STEP 2000010
#define  THERMO 5
#define  CORE 12
#define  OUTSCALE 2
#define  PI 3.1415926

#define  MD_STEP 1   // fs
#define  MD_DUMP 2
#define  MD_TEMP 298 // K


const   double    KB=1.3806504E-23;  // J/K
const   double    U_eV2J=1.602E-19;    
const   double    U_ps2s=1.0E-12;
const   double    U_A2m=1.0E-10;
const   double    U_fs2s=1.0E-15;
const 	double 	  U_atm2Pa=1.01325E5;


double flux[MAX_STEP][MAX_DIM+1];
double autof[MAX_STEP][MAX_DIM+1];
double intg[MAX_STEP][MAX_DIM+1];
double power[MAX_STEP][MAX_DIM+1];


int step;

// double in[MAX_STEP];
// double out[MAX_STEP];


int log_2(int s)
{
	int i, ss;
	i=1;
	ss=s;
	while (ss>1)
	{
		ss=ss/2;	i++;
	}
	return i;
}

int exp_2(int s)
{
	int i, ss;
	i=s;
	ss=1;
	while (i>0)
	{
		ss=ss*2;	i--;
	}
	return ss;
}


void TrapACF(int method)
{
    double      convert, V, scale, box[3];
    FILE        *fbox, *f2;
    int         i, dim;

    convert=U_atm2Pa*U_atm2Pa*U_fs2s*U_A2m*U_A2m*U_A2m;  
    fbox=fopen("info.box", "r");
    fscanf(fbox, "%lf%lf%lf", &box[0], &box[1], &box[2]);
    V=box[0]*box[1]*box[2];
    scale = convert/(KB*MD_TEMP)*MD_STEP*MD_DUMP*V;

    printf("scale = %10.5E\n", scale);

    f2=fopen("03.trap.acf.dat","w");

    // method 1
if (method ==1)
{
    for (dim=0;dim<MAX_DIM;dim++)
    {
        intg[0][dim]=scale*autof[0][dim];
        for (i=1;i<step;i++)
        {
            intg[i][dim]=intg[i-1][dim]+scale*autof[i][dim];
        }
    }
}
    // method 2
if (method==2)
{
    for (dim=0;dim<MAX_DIM;dim++)
    {
        intg[dim][0]=scale*autof[dim][0]/2.0;
        for (i=1;i<step;i++)
        {
            intg[i][dim]=intg[i-1][dim]+scale*(autof[i][dim]+autof[i-1][dim])/2.0;
            // printf("dim = %d i=%d %lf \n", dim, i, autof[i][dim]);
        }
    }
}
    // method 3 simposion
if (method==3)
{
    for (dim=0;dim<MAX_DIM;dim++)
    {
        intg[dim][0]=0.0;
        for (i=2;i<step;i=i+2)
        {
            double ss;
            ss=scale*(autof[i][dim]+4*autof[i-1][dim]+autof[i-2][dim])/3.0;
            intg[i/2][dim]=intg[i/2-1][dim]+ss;
            // printf("dim = %d i=%d %lf \n", dim, i, autof[i][dim]);
        }
    }
}
    
    for (i=0;i<step;i++)
    {
        intg[i][3]=(intg[i][0]+intg[i][1]+intg[i][2])/3.0;
    }

    // fprintf(f2, "nstep  v_ave   v1  v2   v3\n");
    fprintf(f2, "#### t[ps] ### trap {<Q(0)*Q(t)>dt ave xy xz yz\n");

if (method==1 || method==2)
{
    for (i=0;i<step;i++)
    {
        fprintf(f2, "%12.6lf ", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "%20.6E %20.6E %20.6E %20.6E \n", intg[i][3],intg[i][0],intg[i][1],intg[i][2]);  
        // fprintf(f2, "%20.6E\n",intg[i][0]);   
    }
}
if (method==3)
{
    for (i=2;i<step;i=i+2)
    {
        fprintf(f2, "%12.6lf ", i*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "%20.6E %20.6E %20.6E %20.6E \n", intg[i/2][3],intg[i/2][0],intg[i/2][1],intg[i/2][2]);  
        // fprintf(f2, "%20.6E\n",intg[i][0]);   
    }
}

    fclose(f2);    
}


int main()
{
    int i, dim;
    double x,y,z, t;
    double pxx, pyy, pzz, pxy, pyz, pxz;
    double v;
    char header_char[300];
    FILE *f0, *f2;

//  for FFTW
    fftw_complex *in, *out;
    fftw_plan p;

//  open input files 
    f0=fopen("./info.press.data", "r");
    //f0=fopen("../test.data", "r");

    step=0;
    fgets(header_char, 300, f0);
    fgets(header_char, 300, f0);

    while (!feof(f0))
    {
    	step++;
        if (i > 100) exit(0);
     	fscanf(f0, "%lf", &t);
        fscanf(f0, "%lf%lf%lf%lf%lf%lf\n", &pxy, &pxz, &pyz, &pxx, &pyy, &pzz);
    	flux[step-1][0]=pxy;  
        flux[step-1][1]=pxz;  
        flux[step-1][2]=pyz;       
    }
    fclose(f0);

    printf("Total step = %d\n", step);
    int l;
    l=exp_2(log_2(step*2-1));
    printf("Total length = %d\n", l);

    //l=step;	
    for (i=step;i<l;i++)
    {
        flux[step][0]=0;
        flux[step][1]=0;
        flux[step][2]=0;          
    }

    printf("Start FFT\n");
    for (dim=0;dim<MAX_DIM;dim++)
    {
        // FFT
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
        
        for (i=0;i<l;i++)
        {
            in[i][0]=flux[i][dim];
            // printf("input %d %lf\n", i+1,in[i][0]);
            in[i][1]=0.0;
        }
        p = fftw_plan_dft_1d(l, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);
        
        // Inverse FFT
        // printf("Forward FFT\n");
        for (i=0;i<l;i++)
        {
            v=out[i][0]*out[i][0]+out[i][1]*out[i][1];
            in[i][0]=v;
            in[i][1]=0.0;
            //printf("%lf + %lfi  v=%lf\n", out[i][0], out[i][1], v);
        }

        printf("Backward FFT\n");
        p = fftw_plan_dft_1d(l, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);

        for (i=0;i<l;i++)
        {
            autof[i][dim]=out[i][0]/l/(step-i);
            //printf("%E   %E\n", out[i][0]/l, out[i][1]/l);
        }

        fftw_free(in); 
        fftw_free(out);

    }

    // output ACF 
    printf("Output ACF\n");
    double aveacf;
    double scaleacf;

    scaleacf=U_atm2Pa*U_atm2Pa; // unit J/s
    f2=fopen("./0facf.data", "w");
    fprintf(f2, "#### t[ps] ### ACF <P(0)*P(t)> [Pa*Pa]  xy, xz, yz, ave\n");
    printf("%lf\n", MD_STEP*MD_DUMP*1.0E-3);
    for (i=0;i<step;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6E  ", autof[i][0]*scaleacf);
        fprintf(f2, "  %20.6E  ", autof[i][1]*scaleacf);
        fprintf(f2, "  %20.6E  ", autof[i][2]*scaleacf);
        aveacf=(autof[i][0]+autof[i][1]+autof[i][2])/3.0;
        fprintf(f2, "  %20.6E\n", aveacf*scaleacf);
        // fprintf(f2, "  %20.6lf  %20.6lf  %20.6lf ", autof[i][0], autof[i][1], autof[i][2]);
        // aveacf=(autof[i][0]+autof[i][1]+autof[i][2])/3.0;
        // fprintf(f2, "%20.6lf\n", aveacf);
    }
    fclose(f2);

    // output normalized ACF 
    f2=fopen("./01_norm_facf.data", "w");
    fprintf(f2, "#### t[ps] ### ACF <Q(0)*Q(t)> \n");
    for (i=0;i<step;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6lf  ", autof[i][0]/autof[0][0]);
        fprintf(f2, "  %20.6lf  ", autof[i][1]/autof[0][1]);
        fprintf(f2, "  %20.6lf  ", autof[i][2]/autof[0][2]);
        aveacf=(autof[i][0]+autof[i][1]+autof[i][2])/(autof[0][0]+autof[0][1]+autof[0][2])/3.0;
        fprintf(f2, "  %20.6lf\n", aveacf);
    }
    fclose(f2);
    
    // integeral the ACF
    printf("Output integeral of ACF\n");
    printf("Viscosity calculation\n");
    TrapACF(2);
    return 0;
}
