// calculate auto-correlation Function by FFT
// For viscosity calculation 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fftw3.h"


#define  MAX_DIM 3
#define  MAX_STEP 530010
#define  OUTSCALE 10
#define  PI 3.1415926

#define  MD_STEP 1.0   // fs
#define  MD_DUMP 10
#define  MD_TEMP 298 // K

#define MAX_NATOMS 2000


const   double    KB=1.3806504E-23;  // J/K
const   double    U_eV2J=1.602E-19;    
const   double    U_ps2s=1.0E-12;
const   double    U_A2m=1.0E-10;
const   double    U_fs2s=1.0E-15;
const 	double 	  U_atm2Pa=1.01325E5;


short vel[MAX_NATOMS][MAX_STEP]; 
    // velocity of particle
double autof[MAX_STEP];
double sacf[MAX_STEP];


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


int main()
{
    int i, j, k, dim;
    int natoms;
    double x,y,z, t;
    int id, mol, type;
    double  q, x1, y1, z1, xu, yu, zu, vx, vy, vz; 
    double v;
    char header_char[300];
    FILE *f0, *f2;

//  for FFTW
    fftw_complex *in, *out;
    fftw_plan p;

//  open input files 
    f0=fopen("../RCF.info.dat", "r");

    step=0;
    natoms=1580;

    while (!feof(f0))
    {
    	step++;
	if (step>1999900) break;
        for (i=0;i<1;i++) fgets(header_char, 300, f0);
        for (i=0;i<natoms;i++)
        {
            // ITEM: ATOMS id mol type xu yu zu vx vy vz
            fscanf(f0, "%d%d",&id, &vel[i][step-1]);
            fscanf(f0, "\n");         
            // printf("%d\n", vel[i][step-1]);      
        }
    }
    fclose(f0);

    printf("Total step = %d\n", step);
    int l;
    l=exp_2(log_2(step*2-1));
    printf("Total length = %d\n", l);

    //l=step;	
    for (i=0;i<natoms;i++)
    {
        for (j=step;j<l;j++)
        {
            vel[i][j]=0;
            vel[i][j]=0;
            vel[i][j]=0;
        }        
    }

    for (i=0;i<l;i++)
    {
        sacf[i]=0.0;
    }

    printf("Start FFT\n");
    for (k=0;k<natoms;k++)
    {
        // FFT
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
        //printf("Forward FFT %d \n", k);
        for (i=0;i<l;i++)
        {
            in[i][0]=vel[k][i];
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
            // printf("%lf + %lfi  v=%lf\n", out[i][0], out[i][1], v);
        }

        //printf("Backward FFT %d \n", k);
        p = fftw_plan_dft_1d(l, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);

        for (i=0;i<l;i++)
        {
            autof[i]=out[i][0]/l/(step-i);
            // printf("%E   %E\n", out[i][0]/l, out[i][1]/l);
        }
        fftw_free(in); 
        fftw_free(out);
        for (i=0;i<l;i++)
        {
            sacf[i]+=autof[i];
        }
    }

    for (i=0;i<l;i++)
    {
        // sacf[i]/=natoms;
    }


    // output ACF 
    printf("Output ACF\n");
    double scaleacf;

    scaleacf=1.0; // unit m/s * m/2
    f2=fopen("./0rtcf.data", "w");
    fprintf(f2, "#### t[ps] ### RTCF <V(0)*V(t)> \n");
    for (i=0;i<step/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6E  ", sacf[i]*scaleacf);
        fprintf(f2, "\n");
    }
    fclose(f2);

    // output normalized ACF 
    f2=fopen("./01_norm_rtcf.data", "w");
    fprintf(f2, "#### t[ps] ### RTCF <Q(0)*Q(t)> \n");
    for (i=0;i<step/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6E  ", sacf[i]/sacf[0]);
        fprintf(f2, "\n");
    }
    fclose(f2);

    return 0;
}
