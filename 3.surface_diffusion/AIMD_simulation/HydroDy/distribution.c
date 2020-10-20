/// Functional tools 
/// By Ke Zhou 
// ###########################################
// ###########################################
// ###########################################

#include "math_zhou.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define MAX_ATOMS 500
#define PI 3.1415926

double CHAIR2MASS_disp(char element)
{
    if (element=='C') return 12.00;
    if (element=='H') return 1.00;
    if (element=='O') return 16.00;   
    return -1;     
}

int CHAIR2ID_disp(char element)
{
    if (element=='C') return 12;
    if (element=='H') return 1;
    if (element=='O') return 8;   
    return -1;     
}

void HistogramIt_1D(double (*his)[3], double samples, double st, \
    double dt, double fac, int ind)
{
    int kk; 
    kk=(samples-st)/dt;
    his[kk][ind]=his[kk][ind]+1.0*fac;
}
    

void distributions_density(int nwalls, \
    double *box, char *filename)
{
    int    DIM;
    double r[MAX_ATOMS][3];
    double his[800][3];
    // ... 0: density of water 
    // ... 1: Oxygen (O) atom 
    // ... 2: Hydrogen (H) atom 
    char   char_tp;

    int    bins;
    int    len, total_len, i, j, natoms;

    int    tp[MAX_ATOMS];
    double rst, rm, delt, mass[MAX_ATOMS], scale;
    double s[3];
    int    discard_frames;
    
    char header_char[300];
    char output_file[300];
    double data[8];

    FILE *fin;
    FILE *fout;

    DIM=3;
    discard_frames=2000;
    printf(" Calculate density distribution functions: \n");

    rm=20.0;
    rst=0.0;
    bins=400;
    delt=(rm-rst)/bins;
    
    for (i=0;i<bins;i++)
    {
        for (j=0;j<3;j++)
        {
            his[i][j]=0.0; 
        }
    }

    len=0;   
    total_len=0;
    fin=fopen(filename, "r");
  

    while (feof(fin)==0)
    {
        len++;
        if (len % 1000==0) printf("Now is %d\n", len);
        fscanf(fin, "%d\n", &natoms);
        fgets(header_char, 300, fin);
        // puts(header_char);
        // printf("Natoms = %d \n", natoms);
        for (i=0;i<natoms;i++)
        {
            fscanf(fin, "  %c%lf%lf%lf\n", &char_tp, \
                &r[i][0], &r[i][1], &r[i][2]);
            tp[i]=CHAIR2ID_disp(char_tp);
            mass[i]=CHAIR2MASS_disp(char_tp);
            // if (tp==12) continue;
            // HistogramIt_1D(his, r[i][2], rst, delt, mass);
            // printf("%c tp=%d, z=%lf, mass=%lf\n",char_tp, tp, r[i][2], mass);      
        }

        if (len<discard_frames) continue;
        for (i=0;i<natoms;i++)
        {
            // if (tp[i]==12) continue;
           
            if (tp[i]==8)
            {
                HistogramIt_1D(his, r[i][2], rst, delt, 1.00,    1);
                HistogramIt_1D(his, r[i][2], rst, delt, mass[i], 0);    
            }
            if (tp[i]==1)
            {
                HistogramIt_1D(his, r[i][2], rst, delt, 1.00,    2);
                HistogramIt_1D(his, r[i][2], rst, delt, mass[i], 0);    
            }
            // printf("%c tp=%d, z=%lf, mass=%lf\n",char_tp, tp, r[i][2], mass);      
        }        
        total_len++;   
        // printf(" check %d, %lf\n", len, his[40]);  
    }

    printf(" >>The length of input file is %d frames \n", len);
    printf(" >>          selected legth is %d frames \n", total_len);
    len=total_len;
    
    for (i=0;i<bins;i++)
    {
       for (j=0;j<3;j++)
       {
            his[i][j]/=len;
       }
    }

    scale=1.0E27/(6.022E23)/delt/box[0]/box[1];
    // scale=scale/2.0; symetry    
    // cal density
    printf("Scaler = % lf \n", scale);
    for (i=0;i<bins;i++)
    {
        his[i][0]*=scale;
    }   
    // noramlization for O and H distribution
    s[1]=0.0;
    s[2]=0.0;
    for (i=0;i<bins;i++)
    {
        s[1]+=his[i][1];
        s[2]+=his[i][2];
    }   
    for (i=0;i<bins;i++)
    {
        his[i][1]/=s[1];
        his[i][2]/=s[2];
    }   


    sprintf(output_file, "density.profile.dat");
    fout=fopen(output_file, "w");  
    printf(" >> Out put density profile, filename: %s\n", output_file);
    for (i=0;i<bins;i++)
    {
        fprintf(fout, "%12.6lf ", rst+(i+1.0)*delt);
        fprintf(fout, "%12.6lf ", his[i][0]);
        fprintf(fout, "%12.6lf ", his[i][1]*10);
        fprintf(fout, "%12.6lf ", his[i][2]*10);
        fprintf(fout, "\n");
    }
    printf(" >> Please plot the figures using gnuplot\n");
    fclose(fout);
    fclose(fin);
   
}

void distributions_pair_distance(int nwalls, \
    double *box, char *filename)
{
   int    DIM;
    double r[MAX_ATOMS][3];
    double his[400][3];
    // ... 0: O-H, intra 
    // ... 1: O-H, inter 
    // ... 2: O-O, inter 
    char   char_tp;

    int    bins;
    int    len, total_len, i, j, natoms;

    int    tp[MAX_ATOMS];
    double rst, rm, delt;
    double s[3], dr;
    int    discard_frames;
    
    char header_char[300];
    char output_file[300];

    FILE *fin;
    FILE *fout;

    DIM=3;
    discard_frames=2000;
    printf(" Calculate atom-pair distribution functions: \n");

    rm=4.0;
    rst=0.0;
    bins=400;
    delt=(rm-rst)/bins;

    
    for (i=0;i<bins;i++)
    {
        for (j=0;j<3;j++)
        {
            his[i][j]=0.0; 
        }
    }

    len=0;   
    total_len=0;
    fin=fopen(filename, "r");
  

    while (feof(fin)==0)
    {
        len++;
        if (len % 1000==0) printf("Now is %d\n", len);
        fscanf(fin, "%d\n", &natoms);
        fgets(header_char, 300, fin);
        // puts(header_char);
        // printf("Natoms = %d \n", natoms);
        for (i=0;i<natoms;i++)
        {
            fscanf(fin, "  %c%lf%lf%lf\n", &char_tp, \
                &r[i][0], &r[i][1], &r[i][2]);
            tp[i]=CHAIR2ID_disp(char_tp);
        }

        if (len<discard_frames) continue;

        for (i=0;i<natoms;i++)
        {
            for (j=i+1;j<natoms;j++)
            {
                if ((tp[i]==8 && tp[j]==1) || (tp[i]==1 && tp[j]==8))
                {
                    dr=Distant_PBC(r[i], r[j], box);
                    if (dr<1.20)
                    {
                        HistogramIt_1D(his, dr, rst, delt, 1.00, 0);
                    }
                    else if (dr<rm && dr>1.20)
                    {
                        HistogramIt_1D(his, dr, rst, delt, 1.00, 1);
                    }
                }
                if (tp[i]==8 && tp[j]==8)
                {
                    dr=Distant_PBC(r[i], r[j], box);
                    if (dr<rm)
                    {
                        HistogramIt_1D(his, dr, rst, delt, 1.00, 2);
                    }

                }

            }
        }
     
        total_len++;   
        // printf(" check %d, %lf\n", len, his[40]);  
    }

    printf(" >>The length of input file is %d frames \n", len);
    printf(" >>          selected legth is %d frames \n", total_len);
    len=total_len;
    
    for (i=0;i<bins;i++)
    {
       for (j=0;j<3;j++)
       {
            his[i][j]/=len;
       }
    }
 
    // noramlization for O and H distribution
    s[0]=0.0;
    s[1]=0.0;
    s[2]=0.0;
    for (i=0;i<bins;i++)
    {
        s[0]+=his[i][0];
        s[1]+=his[i][1];
        s[2]+=his[i][2];
    }   
    for (i=0;i<bins;i++)
    {
        his[i][0]/=s[0];
        his[i][1]/=s[1];
        his[i][2]/=s[2];
    }   


    sprintf(output_file, "density.pair.profile.dat");
    fout=fopen(output_file, "w");  
    printf(" >> Out put pair-density profile, filename: %s\n", output_file);
    for (i=0;i<bins;i++)
    {
        fprintf(fout, "%12.6lf ", rst+(i+1.0)*delt);
        fprintf(fout, "%12.6lf ", his[i][0]);
        fprintf(fout, "%12.6lf ", his[i][1]);
        fprintf(fout, "%12.6lf ", his[i][2]);
        fprintf(fout, "\n");
    }
    printf(" >> Please plot the figures using gnuplot\n");
    fclose(fout);
    fclose(fin);


}

// void distributions_atomic(int natoms, int nwalls, char *filename, int stp)
// {
    

// }

// void IR_spectra(int natoms, int nwalls, char *filename)
// {
    

// }

// void water_angle(int natoms, int nwalls, char *filename)
// {
    

// }
