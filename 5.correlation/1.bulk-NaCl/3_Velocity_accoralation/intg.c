// calculate auto-correlation Function by FFT
// For viscosity calculation 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// #include "fftw3.h"

// #define OUTPUT_SINGLE_FILE
// #define NOT_OUTPUT_MC
// #define NOT_OUTPUT_MM
// #define NOT_OUTPUT_CC



#define  MAXL 5000
#define  MAX_DIM 3



#define  THERMO 5
#define  PI 3.1415926

#define  MD_STEP 1.0   // fs
#define  MD_DUMP 10
#define  MD_TEMP 298 // K

#define MAX_NATOMS 1


double sacf[MAXL+1][MAX_DIM+1];
double intg[MAXL+1][MAX_DIM+1];

const   double    KB=1.3806504E-23;  // J/K
const   double    U_eV2J=1.602E-19;    
const   double    U_ps2s=1.0E-12;
const   double    U_A2m=1.0E-10;
const   double    U_fs2s=1.0E-15;
const 	double 	  U_atm2Pa=1.01325E5;


void Load_file(FILE *fin, int Len)
{
    int i, j;
    double data[9];
    char header_char[300];

    fgets(header_char, 300, fin);
    for (i=0;i<Len;i++)
    {
        for (j=0;j<9;j++) fscanf(fin, "%lf", &data[j]);
        fscanf(fin, "\n");
        sacf[i][0]=data[1];
        sacf[i][1]=data[2];
        sacf[i][2]=data[3];
        sacf[i][3]=data[4];           
    }
}


void TrapACF(int method, FILE *fout, int Len)
{
    double      convert, V, scale, box[3];
    int         i, dim;

    convert=U_fs2s;  
    scale = convert*MD_STEP*MD_DUMP*1.0E9; // unit 10-^9 m2/s

    printf("scale = %10.5E\n", scale);

    // method 1
if (method ==1)
{
    for (dim=0;dim<MAX_DIM;dim++)
    {
        intg[0][dim]=scale*sacf[0][dim];
        for (i=1;i<Len;i++)
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
        for (i=1;i<Len;i++)
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
        for (i=2;i<Len;i=i+2)
        {
            double ss;
            ss=scale*(sacf[i][dim]+4*sacf[i-1][dim]+sacf[i-2][dim])/3.0;
            intg[i/2][dim]=intg[i/2-1][dim]+ss;
            // printf("dim = %d i=%d %lf \n", dim, i, autof[i][dim]);
        }
    }
}
    
    for (i=0;i<Len;i++)
    {
        intg[i][3]=(intg[i][0]+intg[i][1]+intg[i][2])/3.0;
    }

    // fprintf(f2, "nstep  v_ave   v1  v2   v3\n");
    fprintf(fout, "#### t[ps] ### trap {<Q(0)*Q(t)>dt ave x y z\n");

if (method==1 || method==2)
{
    for (i=0;i<Len;i++)
    {
        fprintf(fout, "%12.6lf ", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(fout, "%20.6E %20.6E %20.6E %20.6E \n", intg[i][3],intg[i][0],intg[i][1],intg[i][2]);  
        // fprintf(f2, "%20.6E\n",intg[i][0]);   
    }
}
if (method==3)
{
    for (i=2;i<Len;i=i+2)
    {
        fprintf(fout, "%12.6lf ", i*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(fout, "%20.6E %20.6E %20.6E %20.6E \n", intg[i/2][3],intg[i/2][0],intg[i/2][1],intg[i/2][2]);  
        // fprintf(f2, "%20.6E\n",intg[i][0]);   
    }
}

  
}




int main()
{
    int i, j, k, k1, k2, dim;
    double x,y,z, t;
	int id, itype, imol;
    int len, l, tim;
    int p_case;
    double chg, xu, yu, zu, vx, vy, vz; 
    double v, scaleacf;
    double rat;
    char header_char[300], filename[300];
    int index_i, index_j;


    FILE *fin, *f2;
    FILE *fsingle;

    fin=fopen("./vacf_self_MM.data", "r");
    f2 =fopen("./vacf_self_MM_intg.data", "w");
    Load_file(fin, 1000);
    TrapACF(1, f2, 1000);
    fclose(fin);
    fclose(f2);
    
    fin=fopen("./vacf_self_CC.data", "r");
    f2 =fopen("./vacf_self_CC_intg.data", "w");
    Load_file(fin, 1000);
    TrapACF(1, f2, 1000);
    fclose(fin);
    fclose(f2);

    fin=fopen("./vacf_cross_MC.data", "r");
    f2 =fopen("./vacf_cross_MC_intg.data", "w");
    Load_file(fin, 1000);
    TrapACF(1, f2, 1000);
    fclose(fin);
    fclose(f2);

    fin=fopen("./vacf_cross_MM.data", "r");
    f2 =fopen("./vacf_cross_MM_intg.data", "w");
    Load_file(fin, 1000);
    TrapACF(1, f2, 1000);
    fclose(fin);
    fclose(f2);

    fin=fopen("./vacf_cross_CC.data", "r");
    f2 =fopen("./vacf_cross_CC_intg.data", "w");
    Load_file(fin, 1000);
    TrapACF(1, f2, 1000);
    fclose(fin);
    fclose(f2);

    return 0;
}
