/// About charge analyse
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

void HistogramIt_charge(double *his, double samples, double st, \
    double dt, double fac)
{
    int kk; 
    kk=(samples-st)/dt;
    his[kk]=his[kk]+1.0*fac;
}
    

void distributions_charge(int nwalls, \
    double *box, char *filename, char *filename1)
{
    int    DIM;
    double r[MAX_ATOMS][3];
    double his[800];
    // ... 0: density of charge 

    char   char_tp;

    int    bins;
    int    len, total_len, i, j, natoms, id, kd;

    int    tp[MAX_ATOMS];
    double rst, rm, delt, muq[MAX_ATOMS], scale;
    double s[3], popu, icharge, sq1, sq2, s1, s2;
    int    discard_frames;
    
    char header_char[300];
    char output_file[300];
    double data[8];

    FILE *finq, *finr;
    FILE *fout;

    DIM=3;
    discard_frames=500;
    printf(" ################################################\n");
    printf(" ################################################\n");
    printf(" Calculate charge density distribution functions: \n");

    rm=20.0;
    rst=0.0;
    bins=400;
    delt=(rm-rst)/bins;
    
    for (i=0;i<bins;i++)
    {
        his[i]=0.0; 
    }

    len=0;   
    total_len=0;
    finq=fopen(filename, "r");
    finr=fopen(filename1, "r");
    
    sq1=0.0; // total charge of the top graphene plane 
    sq2=0.0; // total charge of the buttom graphene plane


    while (feof(finq)==0)
    {
        len++;
        if (len % 100==0) printf("Now is %d\n", len);
        // input position 
        fscanf(finr, "%d\n", &natoms);
        fgets(header_char, 300, finr);
        for (i=0;i<natoms;i++)
        {
            fscanf(finr, "  %c%lf%lf%lf\n", &char_tp, \
                &r[i][0], &r[i][1], &r[i][2]);
        }

        // input Mulliken Population
        s1=0.0;
        s2=0.0;
        // for (i=0;i<5;i++) 
        // {
            // fgets(header_char, 300, finq);
            // puts(header_char);
        // }
        fscanf(finq, "\n");
        fgets(header_char, 300, finq);
        fgets(header_char, 300, finq);
        fscanf(finq, "\n");
        fgets(header_char, 300, finq);
        
        for (i=0;i<natoms;i++)
        {
            fscanf(finq, "%d     %c%d%lf%lf\n", &id, &char_tp, &kd, \
                &popu, &icharge);
            tp[i]=kd;
            muq[i]=icharge;
            // printf("%d, kd= %d,  icharge = %lf\n", id, kd,icharge);
        }
        fgets(header_char, 300, finq);
        fscanf(finq, "\n");
        fgets(header_char, 300, finq);

        // for (i=0;i<1;i++) fgets(header_char, 300, finq);
        // fscanf(finq, "\n");
        // printf("#################\n");

        // if (len>2) break;

        if (len<discard_frames) continue;

        for (i=0;i<nwalls/2;i++) s1+=muq[i];
        for (i=nwalls/2;i<nwalls;i++) s2+=muq[i];

        // printf("XXXXXXXXXXXXXXXXX s1=%lf, s2=%lf\n", s1, s2);
        
        for (i=0;i<natoms;i++)
        {
            HistogramIt_charge(his, r[i][2], rst, delt, muq[i]);
        }        
        total_len++;   
        sq1+=s1;
        sq2+=s2;
        // printf(" check %d, %lf\n", len, his[40]);  
    }

    printf(" >>The length of input file is %d frames \n", len);
    printf(" >>          selected legth is %d frames \n", total_len);
    len=total_len;
    
    for (i=0;i<bins;i++)
    {
        his[i]/=len;
    }
    sq1/=len;
    sq2/=len;
    
    printf(" total chatge on top graphene %lf\n", sq1);
    printf(" total chatge on buttom graphene %lf\n", sq2);
    printf("                         average %lf\n", (sq1-sq2)/2.0);
    
    sprintf(output_file, "charge.density.profile.dat");
    fout=fopen(output_file, "w");  
    printf(" >> Out put charge density profile, filename: %s\n", output_file);
    for (i=0;i<bins;i++)
    {
        fprintf(fout, "%12.6lf ", rst+(i+1.0)*delt);
        fprintf(fout, "%12.6lf ", his[i]);
        fprintf(fout, "\n");
    }
    printf(" >> Please plot the figures using gnuplot\n");
    fclose(fout);
    fclose(finr);
    fclose(finq);  

}

