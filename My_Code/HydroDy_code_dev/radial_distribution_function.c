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

#define MAX_ATOMS 300
#define PI 3.1415926


int CHAIR2ID(char element)
{
    if (element=='C') return 12;
    if (element=='H') return 1;
    if (element=='O') return 8;   
    return -1;     
}
void radial_distribution_functions(int nwalls, \
    double *box, char *filename, int obs1, int obs2)
{
    int    DIM;
    double r[MAX_ATOMS][3], r1[MAX_ATOMS][3], r2[MAX_ATOMS][3];
    int    tp[MAX_ATOMS];
    double his[500][2], his0[500][2];
    char   char_tp;

    double rm, rst, delt, rho, dv, r0;
    int    bins, rr;
    int    len, total_len, i, j, natoms;
    int    step, n1, n2;

    int    discard_frames;
    double rr1, rr2;

    char header_char[300];
    char output_file[300];
    double data[8];

    FILE *fin;
    FILE *fout;

    DIM=3;
    discard_frames=2000;
    printf(" Calculate radial distribution functions g(r): \n");

    rm=8.0;
    rst=0.0;
    bins=300;
    delt=(rm-rst)/bins;
    
    for (i=0;i<bins;i++)
    {
        his[i][0]=0.0;
        his[i][1]=0.0;  
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
        n1=0;
        n2=0;
        // printf("Natoms = %d \n", natoms);
        for (i=0;i<natoms;i++)
        {
            fscanf(fin, "  %c%lf%lf%lf\n", &char_tp, \
                &r[i][0], &r[i][1], &r[i][2]);
            tp[i]=CHAIR2ID(char_tp);
            // printf("%c tp=%d\n",char_tp, tp[i]);
            if (tp[i]==obs1)
            {
                r1[n1][0]=r[i][0];
                r1[n1][1]=r[i][1];
                r1[n1][2]=r[i][2];
                n1++;                   
            }
            if (tp[i]==obs2)
            {
                r2[n2][0]=r[i][0];
                r2[n2][1]=r[i][1];
                r2[n2][2]=r[i][2];
                n2++;                   
            }         
        }


        // printf("N1 =%d N2=%d\n", n1 ,n2);
        if (DIM ==3) rho=n2/box[0]/box[1]/box[2];        
        if (DIM ==2) rho=n2/box[0]/box[1];
        // modify the density if needed
        if (DIM ==3) rho=n2/box[0]/box[1]/(10.0-2.0*2);      
        if (DIM ==3) rho=64/12.42/12.42/12.42; 

        for (i=0;i<bins;i++)
        {
            his0[i][0]=0.0;
            his0[i][1]=0.0; 
        }

        for (i=0;i<n1;i++)
        {
            for (j=0;j<n2;j++)
            {
                if (obs1==obs2 && i==j) continue;
                r0=Distant_PBC(r1[i], r2[j], box)-rst;
                rr=ceil(r0/delt);
                rr1=rr*delt;
                rr2=(rr-1)*delt;
                // if (DIM ==3) dv=4.0/3.0*3.1415926*((r0+delt)*(r0+delt)*(r0+delt)-r0*r0*r0);
                // if (DIM ==2) dv=3.1415926*((r0+delt)*(r0+delt)-r0*r0);
                if (DIM ==3) dv=4.0/3.0*PI*(rr1*rr1*rr1-rr2*rr2*rr2);
                if (DIM ==2) dv=PI*(rr1*rr1-rr2*rr2);
                his0[rr-1][0]+=1.0/dv/rho;
                his0[rr-1][1]+=1.0;
                // if (dis[i][j]<1.0) printf("Find ...\n");
            }
        }
        if (len<discard_frames) continue;
        total_len++;
        for (i=0;i<bins;i++)
        {
            his[i][0]+=(his0[i][0]/n1);
            his[i][1]+=(his0[i][1]/n1);
        }       
    }

    printf(" >>The length of input file is %d frames \n", len);
    printf(" >>          selected legth is %d frames \n", total_len);
    len=total_len;
    for (i=0;i<bins;i++)
    {
        his[i][0]/=len;
        his[i][1]/=len;
    }

    for (i=1;i<bins;i++)
    {
        his[i][1]=his[i][1]+his[i-1][1];
    }

    sprintf(output_file, "rdf.%d.%d.dat", obs1, obs2);
    fout=fopen(output_file, "w");  
    printf(" >> Out put RDF, filename: %s\n", output_file);
    for (i=0;i<bins;i++)
    {
        fprintf(fout, "%12.6lf ", rst+(i+1.0)*delt);
        fprintf(fout, "%12.6lf ", his[i][0]);
        fprintf(fout, "%12.6lf ", his[i][1]);
        fprintf(fout, "\n");
    }
    printf(" >> Please plot the figures using gnuplot\n");
    fclose(fout);
    fclose(fin);

}

// void mean_square_distance(int natoms, int nwalls, char *filename)
// {
    

// }



// void distributions_atomic(int natoms, int nwalls, char *filename, int stp)
// {
    

// }

// void IR_spectra(int natoms, int nwalls, char *filename)
// {
    

// }

// void water_angle(int natoms, int nwalls, char *filename)
// {
    

// }