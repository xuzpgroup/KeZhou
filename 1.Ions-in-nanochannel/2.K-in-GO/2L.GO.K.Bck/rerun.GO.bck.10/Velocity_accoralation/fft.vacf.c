// calculate auto-correlation Function by FFT
// For viscosity calculation 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fftw3.h"


#define  MAX_DIM 3
#define  MAX_STEP 2100010
#define  THERMO 5
#define  OUTSCALE 100
#define  PI 3.1415926

#define  MD_STEP 1.0   // fs
#define  MD_DUMP 10
#define  MD_TEMP 298 // K

#define MAX_NATOMS 1


const   double    KB=1.3806504E-23;  // J/K
const   double    U_eV2J=1.602E-19;    
const   double    U_ps2s=1.0E-12;
const   double    U_A2m=1.0E-10;
const   double    U_fs2s=1.0E-15;
const 	double 	  U_atm2Pa=1.01325E5;


double vel[MAX_NATOMS][MAX_STEP][3]; 
    // velocity of particle
double autof[MAX_STEP][MAX_DIM+1];
double sacf[MAX_STEP][MAX_DIM+1];
double intg[MAX_STEP][MAX_DIM+1];


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

    convert=U_fs2s;  
    scale = convert*MD_STEP*MD_DUMP*1.0E9; // unit 10^9 m2/s

    printf("scale = %10.5E\n", scale);

    f2=fopen("03.trap.acf.dat","w");

    // method 1
if (method ==1)
{
    for (dim=0;dim<MAX_DIM;dim++)
    {
        intg[0][dim]=scale*sacf[0][dim];
        for (i=1;i<step;i++)
        {
            intg[i][dim]=intg[i-1][dim]+scale*sacf[i][dim];
        }
    }
}
    // method 2
if (method==2)
{
    for (dim=0;dim<MAX_DIM;dim++)
    {
        intg[dim][0]=scale*sacf[dim][0]/2.0;
        for (i=1;i<step;i++)
        {
            intg[i][dim]=intg[i-1][dim]+scale*(sacf[i][dim]+sacf[i-1][dim])/2.0;
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
            ss=scale*(sacf[i][dim]+4*sacf[i-1][dim]+sacf[i-2][dim])/3.0;
            intg[i/2][dim]=intg[i/2-1][dim]+ss;
            // printf("dim = %d i=%d %lf \n", dim, i, autof[i][dim]);
        }
    }
}
    
    for (i=0;i<step;i++)
    {
        intg[i][3]=(intg[i][0]+intg[i][1])/2.0;
    }

    // fprintf(f2, "nstep  v_ave   v1  v2   v3\n");
    fprintf(f2, "#### t[ps] ### trap {<Q(0)*Q(t)>dt ave x y z\n");

if (method==1 || method==2)
{
    for (i=0;i<step/OUTSCALE;i++)
    {
        fprintf(f2, "%12.6lf ", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "%20.6E %20.6E %20.6E %20.6E \n", intg[i][3],intg[i][0],intg[i][1],intg[i][2]);  
        // fprintf(f2, "%20.6E\n",intg[i][0]);   
    }
}
if (method==3)
{
    for (i=2;i<step/OUTSCALE;i=i+2)
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
    f0=fopen("../ions.dat", "r");

    step=0;

    while (!feof(f0))
    {
    	step++;
        for (i=0;i<3;i++) fgets(header_char, 300, f0);
        fscanf(f0, "%d\n", &natoms);
        for (i=0;i<5;i++) fgets(header_char, 300, f0);

        for (i=0;i<natoms;i++)
        {
            // ITEM: ATOMS id mol type xu yu zu vx vy vz
            fscanf(f0, "%d%d%d",&id, &mol, &type);
            fscanf(f0, "%lf%lf%lf", &xu, &yu, &zu);
            fscanf(f0, "%lf%lf%lf", &vx, &vy, &vz);
            fscanf(f0, "\n");               
            if (type==9)
            {
                vel[i][step-1][0]=vx*1.0E5; // [A/fs] to [m/s]
                vel[i][step-1][1]=vy*1.0E5; //
                vel[i][step-1][2]=vz*1.0E5; //
                // printf("vx=%lf vy=%lf vz=%lf \n", vel[i][step-1][0], vel[i][step-1][1], vel[i][step-1][2]);
            }
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
            vel[i][j][0]=0;
            vel[i][j][1]=0;
            vel[i][j][2]=0;
        }        
    }

    for (i=0;i<l;i++)
    {
        for (j=0;j<MAX_DIM+1;j++)
        {
            sacf[i][j]=0.0;
        }
    }

    printf("Start FFT\n");
    for (k=0;k<natoms;k++)
    {
        for (dim=0;dim<MAX_DIM;dim++)
        {
            // FFT
            in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
            out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
            printf("Forward FFT\n");
            for (i=0;i<l;i++)
            {
                in[i][0]=vel[k][i][dim];
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
            for (i=0;i<l;i++)
            {
                sacf[i][dim]+=autof[i][dim];
            }
        }

    }
    for (i=0;i<l;i++)
    {
        for (j=0;j<MAX_DIM;j++)
        {
            sacf[i][j]/=natoms;
        }
        sacf[i][MAX_DIM]=sacf[i][0]+sacf[i][1];
    }


    // output ACF 
    printf("Output ACF\n");
    double scaleacf;

    scaleacf=1.0; // unit m/s * m/2
    f2=fopen("./0vacf.data", "w");
    fprintf(f2, "#### t[ps] ### VCF <V(0)*V(t)> [(m/s)^2]  Vx, Vy, Vz, Vxy \n");
    printf("%lf\n", MD_STEP*MD_DUMP*1.0E-3);
    for (i=0;i<step/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6E  ", sacf[i][0]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[i][1]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[i][2]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[i][3]*scaleacf);
        fprintf(f2, "\n");
    }
    fclose(f2);

    // output normalized ACF 
    f2=fopen("./01_norm_vacf.data", "w");
    fprintf(f2, "#### t[ps] ### VACF <Q(0)*Q(t)> \n");
    for (i=0;i<step/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6lf  ", sacf[i][0]/sacf[0][0]);
        fprintf(f2, "  %20.6lf  ", sacf[i][1]/sacf[0][1]);
        fprintf(f2, "  %20.6lf  ", sacf[i][2]/sacf[0][2]);
        fprintf(f2, "  %20.6lf  ", sacf[i][3]/sacf[0][3]);
        fprintf(f2, "\n");
    }
    fclose(f2);
    
    // integeral the ACF
    printf("Output integeral of ACF\n");
    printf("Diffusion calculation\n");
    TrapACF(3);
    return 0;
}
