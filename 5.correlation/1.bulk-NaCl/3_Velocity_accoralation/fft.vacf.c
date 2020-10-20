// calculate auto-correlation Function by FFT
// For viscosity calculation 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fftw3.h"

// #define OUTPUT_SINGLE_FILE
// #define NOT_OUTPUT_MC
// #define NOT_OUTPUT_MM
// #define NOT_OUTPUT_CC



#define MAXNATOM 100
// #define MAXT 1048576
#define MAXT  524288  // 2<<19, for 5ns

// const int MAXT = 1<<20;  //  for 10ns

// #define MAXT  524288  // 2<<19, for 5ns

#define OUTSCALE 50

#define  MAX_DIM 3
#define  THERMO 5
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

// double 	vM[MAXNATOM*2][MAXT][3];
double  vM[2][MAXT+1][3];
// double 	vC[MAXNATOM][MAXT*2][3];
    // velocity of particle


double autof[MAXT+1][MAX_DIM+1];
double sacf[5][MAXT+1][MAX_DIM+1];
// ... 0 Na-Na self-correlation
// ... 1 Cl-Cl self-correlation
// ... 2 Na-Cl cross-correlation
// ... 3 Na-Na cross-correlation
// ... 4 Cl-Cl cross-correlation



// double sacf1[MAXT+1][MAX_DIM+1];
// double sacf2[MAXT+1][MAX_DIM+1];
// double intg[MAXT][MAX_DIM+1];

double box[3];
int natoms;

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


// void TrapACF(int method)
// {
//     double      convert, V, scale, box[3];
//     FILE        *fbox, *f2;
//     int         i, dim;

//     convert=U_fs2s;  
//     scale = convert*MD_STEP*MD_DUMP*1.0E9; // unit 10^9 m2/s

//     printf("scale = %10.5E\n", scale);

//     f2=fopen("03.trap.acf.dat","w");

//     // method 1
// if (method ==1)
// {
//     for (dim=0;dim<MAX_DIM;dim++)
//     {
//         intg[0][dim]=scale*sacf[0][dim];
//         for (i=1;i<step;i++)
//         {
//             intg[i][dim]=intg[i-1][dim]+scale*sacf[i][dim];
//         }
//     }
// }
//     // method 2
// if (method==2)
// {
//     for (dim=0;dim<MAX_DIM;dim++)
//     {
//         intg[dim][0]=scale*sacf[dim][0]/2.0;
//         for (i=1;i<step;i++)
//         {
//             intg[i][dim]=intg[i-1][dim]+scale*(sacf[i][dim]+sacf[i-1][dim])/2.0;
//             // printf("dim = %d i=%d %lf \n", dim, i, autof[i][dim]);
//         }
//     }
// }
//     // method 3 simposion
// if (method==3)
// {
//     for (dim=0;dim<MAX_DIM;dim++)
//     {
//         intg[dim][0]=0.0;
//         for (i=2;i<step;i=i+2)
//         {
//             double ss;
//             ss=scale*(sacf[i][dim]+4*sacf[i-1][dim]+sacf[i-2][dim])/3.0;
//             intg[i/2][dim]=intg[i/2-1][dim]+ss;
//             // printf("dim = %d i=%d %lf \n", dim, i, autof[i][dim]);
//         }
//     }
// }
    
//     for (i=0;i<step;i++)
//     {
//         intg[i][3]=(intg[i][0]+intg[i][1])/2.0;
//     }

//     // fprintf(f2, "nstep  v_ave   v1  v2   v3\n");
//     fprintf(f2, "#### t[ps] ### trap {<Q(0)*Q(t)>dt ave x y z\n");

// if (method==1 || method==2)
// {
//     for (i=0;i<step/OUTSCALE;i++)
//     {
//         fprintf(f2, "%12.6lf ", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
//         fprintf(f2, "%20.6E %20.6E %20.6E %20.6E \n", intg[i][3],intg[i][0],intg[i][1],intg[i][2]);  
//         // fprintf(f2, "%20.6E\n",intg[i][0]);   
//     }
// }
// if (method==3)
// {
//     for (i=2;i<step/OUTSCALE;i=i+2)
//     {
//         fprintf(f2, "%12.6lf ", i*MD_STEP*MD_DUMP*1.0E-3);
//         fprintf(f2, "%20.6E %20.6E %20.6E %20.6E \n", intg[i/2][3],intg[i/2][0],intg[i/2][1],intg[i/2][2]);  
//         // fprintf(f2, "%20.6E\n",intg[i][0]);   
//     }
// }

//     fclose(f2);    
// }
void inti_()
{
    int i, j, k;
    for (i=0;i<5;i++)
    {
        for (j=0;j<MAXT+1;j++)
        {
            for (k=0;k<MAX_DIM+1;k++)
            {
                sacf[i][j][k]=0.0;
            }

        }
    }
}

int pair_cases(int p, int q)
{
    if (p==q)
    {
        if (p<natoms/2) return 0;
        return 1;
    }
    if (p<natoms/2 && q<natoms/2)
    {
        return 3;
    }
    if (p>=natoms/2 && q>=natoms/2)
    {
        return 4;
    }
    return 2;
}

void inp_basic()
{
    FILE *fin;
    int i;
    double  xlo, xhi, ylo, yhi, zlo, zhi;
    char header_char[300];

    fin=fopen("../ions.dat", "r");

    for (i=0;i<3;i++) fgets(header_char, 300, fin);    
    fscanf(fin, "%d\n", &natoms);
    fgets(header_char, 300, fin);    
    
    fscanf(fin, "%lf%lf\n", &xlo, &xhi);
    fscanf(fin, "%lf%lf\n", &ylo, &yhi);
    fscanf(fin, "%lf%lf\n", &zlo, &zhi);
    fgets(header_char, 300, fin); 
    
    box[0]=xhi-xlo;
    box[1]=yhi-ylo;
    box[2]=zhi-zlo;
    printf(" N_Na = %d, N_Cl =%d. \n", natoms/2, natoms/2);
    printf(" Box_x= %lf \n", box[0]);
    printf(" Box_y= %lf \n", box[1]);
    printf(" Box_z= %lf \n", box[2]);  
    printf("############################################\n"); 

    fclose(fin);
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

//  for FFTW
    fftw_complex *in, *out, *in1, *out1;
    fftw_plan p, p1;

//  open input files 

    printf(" Input basic information ...\n");
    inp_basic();

    printf(" Input velocity ...\n");

    inti_();
    tim=0;
    // return 0;
    printf("Start .....\n");
   
for (index_i=0;index_i<natoms;index_i++)
    for (index_j=index_i;index_j<natoms;index_j++)
    {
    	tim++;
        p_case=pair_cases(index_i, index_j);
////////////////////////////////
#ifdef  NOT_OUTPUT_MC
        if ( p_case==2) continue;
#endif 
#ifdef  NOT_OUTPUT_MM
        if ( p_case==3) continue;
#endif  
#ifdef  NOT_OUTPUT_CC
        if ( p_case==4) continue;
#endif  
        
////////////////////////////////
    fin=fopen("../ions.dat", "r");
    len=0;

    printf("### %d Input velocity %d x %d\n", tim, index_i+1, index_j+1);
    rat=tim/((natoms*(natoms+1))/2+0.000001);
    printf("### ratio = %lf \n", rat);
    printf("[");
    for (i=0;i<floor(64*rat);i++) printf("*");
    for (i=0;i<floor(64-64*rat);i++) printf("_");
    printf("]\n");
    while (!feof(fin))
    {
        if (len>MAXT+1) break;   
        len++;
        if (len % 100000 ==1) printf("Now is %d, step = %d\n", len, step);
    
        // input water and ions from FILE: ions.dat
        for (i=0;i<9;i++) fgets(header_char, 300, fin);

        for (i=0;i<natoms;i++)
        {
            // ITEM: ATOMS id type xu yu zu
            fscanf(fin, "%d%d%d%lf", &id, &imol, &itype, &chg);
            fscanf(fin, "%lf%lf%lf", &vx, &vy, &vz);
            fscanf(fin, "%lf%lf%lf\n", &xu, &yu, &zu);
            if (i==index_i)
            {
                vM[0][len-1][0]=vx*1.0E5; // [A/fs] to [m/s];
                vM[0][len-1][1]=vy*1.0E5; // [A/fs] to [m/s];
                vM[0][len-1][2]=vz*1.0E5; // [A/fs] to [m/s];
            }    
            if (i==index_j)
            {
                vM[1][len-1][0]=vx*1.0E5; // [A/fs] to [m/s];
                vM[1][len-1][1]=vy*1.0E5; // [A/fs] to [m/s];
                vM[1][len-1][2]=vz*1.0E5; // [A/fs] to [m/s];
            }                    
        }
    }
    fclose(fin);
    // printf("### Input velocity %d x %d\n", index_i+1, index_j+1);
    printf("### the length of data %d \n", len);

////////////////////////////////
    printf("### Calculate correlation for case %d\n",  p_case);

    l=MAXT;
    printf("### Total length = %d\n", l);
    step=len;   

    for (i=0;i<2;i++)
    {
        for (j=step;j<l;j++)
        {
            vM[i][j][0]=0;
            vM[i][j][1]=0;
            vM[i][j][2]=0;
        }        
    }

    printf("### Start FFT \n");
    for (dim=0;dim<MAX_DIM;dim++)
    {
        // FFT for M
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
        // printf("Forward FFT %d-th for ions\n", k+1);
        for (i=0;i<l;i++)
        {
            in[i][0]=vM[0][i][dim];
            // printf("input %d %lf\n", i+1,in[i][0]);
            in[i][1]=0.0;
        }
        p = fftw_plan_dft_1d(l, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);

        // FFT for C
        in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
        out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
        // printf("Forward FFT %d-th for ions\n", k+1);
        for (i=0;i<l;i++)
        {
            in1[i][0]=vM[1][i][dim];
            // printf("input %d %lf\n", i+1,in[i][0]);
            in1[i][1]=0.0;
        }
        p1 = fftw_plan_dft_1d(l, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p1);
        fftw_destroy_plan(p1);

        // Inverse FFT
        // printf("Forward FFT\n");
        for (i=0;i<l;i++)
        {
            v=out[i][0]*out1[i][0]+out[i][1]*out1[i][1];
            in[i][0]=v;
            in[i][1]=0.0;
            // printf("%lf + %lfi  v=%lf\n", out[i][0], out[i][1], v);
        }

        // printf("Backward FFT\n");
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
        fftw_free(in1); 
        fftw_free(out1);
        for (i=0;i<l;i++)
        {
            sacf[p_case][i][dim]+=autof[i][dim];
        }
    }
    printf("### End FFT \n\n\n");

    }




// #########################################
// #########################################
// #########################################
// #########################################
// #########################################
// #########################################
    // output ACF 
    printf("\n Output ACF\n");


    for (i=0;i<l;i++)
    {
        for (j=0;j<MAX_DIM;j++)
        {
            for (k=0;k<5;k++)
            {
                sacf[k][i][j]/=(natoms/2);
            } 
        }
        
        for (k=0;k<5;k++)
        {
            sacf[k][i][j]/=(natoms/2);
            sacf[k][i][MAX_DIM]=sacf[k][i][0]+sacf[k][i][1]+sacf[k][i][2];
        }         
    }




    // output ACF 
    printf("\n Output Cross ACF\n");

    scaleacf=1.0; // unit (m/s)^2

    f2=fopen("./vacf_self_MM.data", "w");
    fprintf(f2, "#### t[ps] ### VCF <V(0)*V(t)> [(m/s)^2]  Vx, Vy, Vz, Vxyz \n");
    printf("%lf\n", MD_STEP*MD_DUMP*1.0E-3);
    for (i=0;i<len/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6E  ", sacf[0][i][0]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[0][i][1]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[0][i][2]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[0][i][3]*scaleacf);
        // for normalization
        fprintf(f2, "  %20.6E  ", sacf[0][i][0]/sacf[0][0][0]);
        fprintf(f2, "  %20.6E  ", sacf[0][i][1]/sacf[0][0][1]);
        fprintf(f2, "  %20.6E  ", sacf[0][i][2]/sacf[0][0][2]);
        fprintf(f2, "  %20.6E  ", sacf[0][i][3]/sacf[0][0][3]);
        fprintf(f2, "\n");
    }
    fclose(f2);

    f2=fopen("./vacf_self_CC.data", "w");
    fprintf(f2, "#### t[ps] ### VCF <V(0)*V(t)> [(m/s)^2]  Vx, Vy, Vz, Vxyz \n");
    printf("%lf\n", MD_STEP*MD_DUMP*1.0E-3);
    for (i=0;i<len/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6E  ", sacf[1][i][0]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[1][i][1]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[1][i][2]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[1][i][3]*scaleacf);
        // for normalization
        fprintf(f2, "  %20.6E  ", sacf[1][i][0]/sacf[1][0][0]);
        fprintf(f2, "  %20.6E  ", sacf[1][i][1]/sacf[1][0][1]);
        fprintf(f2, "  %20.6E  ", sacf[1][i][2]/sacf[1][0][2]);
        fprintf(f2, "  %20.6E  ", sacf[1][i][3]/sacf[1][0][3]);
        fprintf(f2, "\n");
    }
    fclose(f2);


    f2=fopen("./vacf_cross_MC.data", "w");
    fprintf(f2, "#### t[ps] ### VCF <V(0)*V(t)> [(m/s)^2]  Vx, Vy, Vz, Vxyz \n");
    printf("%lf\n", MD_STEP*MD_DUMP*1.0E-3);
    for (i=0;i<len/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6E  ", sacf[2][i][0]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[2][i][1]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[2][i][2]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[2][i][3]*scaleacf);
        // for normalization
        fprintf(f2, "  %20.6E  ", sacf[2][i][0]/sacf[2][0][0]);
        fprintf(f2, "  %20.6E  ", sacf[2][i][1]/sacf[2][0][1]);
        fprintf(f2, "  %20.6E  ", sacf[2][i][2]/sacf[2][0][2]);
        fprintf(f2, "  %20.6E  ", sacf[2][i][3]/sacf[2][0][3]);
        fprintf(f2, "\n");
    }
    fclose(f2);

    f2=fopen("./vacf_cross_MM.data", "w");
    fprintf(f2, "#### t[ps] ### VCF <V(0)*V(t)> [(m/s)^2]  Vx, Vy, Vz, Vxyz \n");
    printf("%lf\n", MD_STEP*MD_DUMP*1.0E-3);
    for (i=0;i<len/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6E  ", sacf[3][i][0]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[3][i][1]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[3][i][2]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[3][i][3]*scaleacf);
        // for normalization
        fprintf(f2, "  %20.6E  ", sacf[3][i][0]/sacf[3][0][0]);
        fprintf(f2, "  %20.6E  ", sacf[3][i][1]/sacf[3][0][1]);
        fprintf(f2, "  %20.6E  ", sacf[3][i][2]/sacf[3][0][2]);
        fprintf(f2, "  %20.6E  ", sacf[3][i][3]/sacf[3][0][3]);
        fprintf(f2, "\n");
    }
    fclose(f2);
    

    f2=fopen("./vacf_cross_CC.data", "w");
    fprintf(f2, "#### t[ps] ### VCF <V(0)*V(t)> [(m/s)^2]  Vx, Vy, Vz, Vxyz \n");
    printf("%lf\n", MD_STEP*MD_DUMP*1.0E-3);
    for (i=0;i<len/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+1)*MD_STEP*MD_DUMP*1.0E-3);
        fprintf(f2, "  %20.6E  ", sacf[4][i][0]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[4][i][1]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[4][i][2]*scaleacf);
        fprintf(f2, "  %20.6E  ", sacf[4][i][3]*scaleacf);
        // for normalization
        fprintf(f2, "  %20.6E  ", sacf[4][i][0]/sacf[4][0][0]);
        fprintf(f2, "  %20.6E  ", sacf[4][i][1]/sacf[4][0][1]);
        fprintf(f2, "  %20.6E  ", sacf[4][i][2]/sacf[4][0][2]);
        fprintf(f2, "  %20.6E  ", sacf[4][i][3]/sacf[4][0][3]);
        fprintf(f2, "\n");
    }
    fclose(f2);
    return 0;
}
