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
#define SQRT3 1.732050808

#define ROO 3.60
#define ROH 2.45

void HistogramIt_ang(double (*his)[2], double samples, double st, \
    double dt, double fac, int ind)
{
    int kk; 
    kk=(samples-st)/dt;
    his[kk][ind]=his[kk][ind]+1.0*fac;
}
    
void HistogramIt_ang_2D(double (*his)[800], \
    double samples1, double st1, double dt1, \
    double samples2, double st2, double dt2, double fac)
{
    int kk1, kk2; 
    kk1=(samples1-st1)/dt1;
    kk2=(samples2-st2)/dt2;
    his[kk1][kk2]=his[kk1][kk2]+1.0*fac;
}

void HB_Rotation(int nwalls, \
    double *box, char *filename)
{
    int    DIM;
    double r[MAX_ATOMS][3];
    double his_hb_ang[800][2];
    // ... 0: ang to bottom 
    // ... 1: ang to top 
    
    double his_2D[800][800];
    double s[2], ss;
    
    char   char_tp;
    int    len, total_len, i, j, natoms;

    int    tp[MAX_ATOMS];
    double mass[MAX_ATOMS], scale;
    double rst, rm, delt, rstz, rmz, deltz;
    int    bins, binsz;
    int    H1, H2, O1;
    double vec_OH_1[3], vec_OH_2[3], theta1, theta2;
    int    discard_frames;

    double dHB, aHB, dOO, dOH1, dOH2, sNhB;
    double dir2[3], dir1[3], cosang;
    int    idO, idH1, idH2, jdO, jdH1, jdH2;
    int    flag;
    
    char header_char[300];
    char output_file[300];
    double data[8];

    FILE *fin;
    FILE *fout;

    DIM=3;
    discard_frames=2000;
    printf(" Calculate hydrogen bond distribution: \n");

    rm=180.0;
    rst=0.0;
    bins=100;
    delt=(rm-rst)/bins;
    
    rmz=10.0;
    rstz=0.0;
    binsz=100;
    deltz=(rmz-rstz)/binsz;
    
    
    for (i=0;i<bins;i++)
    {
        for (j=0;j<2;j++)
        {
            his_hb_ang[i][j]=0.0; 
        }
    }
    for (i=0;i<bins;i++)
    {
        for (j=0;j<binsz;j++)
        {
            his_2D[i][j]=0.0;
        }
    }

    len=0;   
    total_len=0;
    fin=fopen(filename, "r");
  
    sNhB=0.0;

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
            tp[i]=CHAIR2ID_math(char_tp);
            mass[i]=CHAIR2MASS_math(char_tp);
            // if (tp==12) continue;
            // HistogramIt_1D(his, r[i][2], rst, delt, mass);
            // printf("%c tp=%d, z=%lf, mass=%lf\n",char_tp, tp, r[i][2], mass);      
        }

        if (len<discard_frames) continue;

        // calculate O-H angle distribution
        for (i=0;i<natoms;i=i+3)
        {
            H1=i;
            H2=i+1;
            O1=i+2;
            Vec_Minus_nromal(r[H1], r[O1], vec_OH_1);
            Vec_Minus_nromal(r[H2], r[O1], vec_OH_2);
            theta1=acos(vec_OH_1[2])*180.0/3.1415926;
            theta2=acos(vec_OH_2[2])*180.0/3.1415926;
         
            if (r[O1][2]<5.0)
            {
                HistogramIt_ang(his_hb_ang, theta1, rst,delt, 1.00, 0);
                HistogramIt_ang(his_hb_ang, theta2, rst,delt, 1.00, 0);
            }
            else 
            {
                HistogramIt_ang(his_hb_ang, theta1, rst,delt, 1.00, 1);
                HistogramIt_ang(his_hb_ang, theta2, rst,delt, 1.00, 1);
            }
            HistogramIt_ang_2D(his_2D, 
                theta1, rst, delt, \
                r[O1][2], rstz, deltz, 1.00);
            HistogramIt_ang_2D(his_2D, 
                theta2, rst, delt, \
                r[O1][2], rstz, deltz, 1.00);
        }        
        total_len++;  

        // calculate H bond 
        // Criterion for H bond 
        // ROO<3.60 ROH<2.45 for SPC/E
        // <H--O..O < 30.0
        // #########################
        dHB=0.0;
        for (i=0;i<natoms/3;i++)
        {
            aHB=0.0;
            // center mol i
            idH1=i*3+0;
            idH2=i*3+1;
            idO =i*3+2;
            // find donor 
            for (j=0;j<natoms/3;j++)
            if (i!=j)
            {
                jdO =j*3+2;
                jdH1=j*3+0;
                jdH2=j*3+1;             
                dOO=Distant_PBC(r[idO], r[jdO], box);
                // printf("dOO %lf\n", dOO);
                if (dOO>ROO) continue;
                dOH1=Distant_PBC(r[idO], r[jdH1], box);
                dOH2=Distant_PBC(r[idO], r[jdH2], box);
                flag=0;
                if (dOH1<ROH) flag=1;
                if (dOH2<ROH) flag=2;
                // printf("dOH %lf %lf\n", dOH1, dOH2);
                if (dOH1>ROH && dOH2>ROH) continue;
                // Vec_Minus_nromal(r[idO], r[jdO], dir1);
                Vec_Minus_PBC(r[idO], r[jdO], dir1, 1, box);
                if (flag==1) 
                {
                    // Vec_Minus_nromal(r[jdH1], r[jdO], dir2);
                    Vec_Minus_PBC(r[jdH1], r[jdO], dir2, 1, box);
                    // printf("For ang\n");
                    cosang=DotProduct(dir1, dir2);
                    if (cosang>SQRT3/2.0)
                    {
                        sNhB+=1.0;
                        dHB+=1.0;
                        aHB+=1.0;
                        // printf("get it %d %d\n", i,j);
                        // histogram 
                    }           
                    continue;   
                }
                if (flag==2)
                {
                    // Vec_Minus_nromal(r[jdH2], r[jdO], dir2);
                    Vec_Minus_PBC(r[jdH2], r[jdO], dir2, 1, box);
                    cosang=DotProduct(dir1, dir2);
                    if (cosang>SQRT3/2.0)
                    {
                        sNhB+=1.0;
                        dHB+=1.0;
                        aHB+=1.0;
                        // printf("get it %d %d\n", i,j);
                        // histogram
                    }               
                    continue;   
                }

            }
            // printf("MOL %d donor %lf\n", i+1, dHB);
            // find acceptor 
            for (j=0;j<natoms/3;j++)
            if (i!=j)
            {
                jdO=j*3+2;      
                dOO=Distant_PBC(r[idO], r[jdO], box);
                if (dOO>ROO) continue;
                dOH1=Distant_PBC(r[jdO], r[idH1], box);
                dOH2=Distant_PBC(r[jdO], r[idH2], box);
                flag=0;
                if (dOH1<ROH) flag=1;
                if (dOH2<ROH) flag=2;
                if (dOH1>ROH && dOH2>ROH) continue;
                // printf("For ang\n");
                Vec_Minus_PBC(r[jdO], r[idO], dir1, 1, box);
                if (flag==1)
                {
                    Vec_Minus_PBC(r[idH1], r[idO], dir2, 1, box);
                    cosang=DotProduct(dir1, dir2);
                    if (cosang>SQRT3/2.0)
                    {
                        sNhB+=1.0;
                        dHB+=1.0;
                        aHB+=1.0;
                        // histogram 
                    }               
                    continue;   
                }
                if (flag==2)
                {
                    Vec_Minus_PBC(r[idH2], r[idO], dir2, 1, box);
                    cosang=DotProduct(dir1, dir2);
                    if (cosang>SQRT3/2.0)
                    {
                        sNhB+=1.0;
                        dHB+=1.0;
                        aHB+=1.0;
                        // histogram
                    }               
                    continue;   
                }

            }
            // printf("acceptor %lf\n", dHB);
        }
        // printf("%d HB %lf  nmol=%lf\n", total_len, 1.0*dHB/natoms*3.0, natoms/3.0);
        if (len % 1000 ==0) printf("Now is %dth trajectory\n", total_len);
        if (len % 1000 ==0) printf("  N water = %d\n", natoms /3);
        if (len % 1000 ==0) printf("  NhB = %lf\n", sNhB*1.0/(natoms/3)/total_len);
        // #########################




        // printf(" check %d, %lf\n", len, his[40]);  
    }
    fclose(fin);

    printf(" >>The length of input file is %d frames \n", len);
    printf(" >>          selected legth is %d frames \n", total_len);
    len=total_len;
    
 
    s[0]=0.0;
    s[1]=0.0;
    ss=0.0;

    for (i=0;i<bins;i++)
    {
       for (j=0;j<2;j++)
       {
            his_hb_ang[i][j]/=len;
            s[j]+=his_hb_ang[i][j];
       }
    }
    for (i=0;i<bins;i++)
    {
        for (j=0;j<binsz;j++)
        {
            his_2D[i][j]/=len;
            ss+=his_2D[i][j];
        }
    }   

    // scale=scale/2.0; symetry    
    // cal density

    // noramlization for O-H distribution
 
    for (i=0;i<bins;i++)
    {
        his_hb_ang[i][0]/=s[0];
        his_hb_ang[i][1]/=s[1];
    }   
    for (i=0;i<bins;i++)
    {
        for (j=0;j<binsz;j++)
        {
            his_2D[i][j]/=ss;
        }
    }   

    sprintf(output_file, "O-H.ang.dat");
    fout=fopen(output_file, "w");  
    printf(" >> Out put O-H angle distribution profile, filename: %s\n", output_file);
    
    fprintf(fout, "%12.6lf %12.6lf %12.6lf \n", rst, 0.0, 0.0);
    for (i=0;i<bins;i++)
    {
        fprintf(fout, "%12.6lf ", rst+(i+1.0)*delt);
        fprintf(fout, "%12.6lf ", his_hb_ang[i][0]);
        fprintf(fout, "%12.6lf ", his_hb_ang[i][1]);
        fprintf(fout, "\n");
    }
    fprintf(fout, "%12.6lf %12.6lf %12.6lf \n", rm+delt, 0.0, 0.0);
    fclose(fout);

    sprintf(output_file, "O-H.ang.z.dat");
    fout=fopen(output_file, "w");      
    // for (j=0;j<binsz;j++)
    // {
    //     fprintf(fout, "%12.6lf ", rst+(0+0.0)*delt);
    //     fprintf(fout, "%12.6lf ", rstz+(j+1.0)*deltz);
    //     fprintf(fout, "%12.10lf ", 0.00);
    //     fprintf(fout, "\n");     	
    // }
    fprintf(fout, "\n"); 
    for (i=0;i<bins;i++)
    {
        for (j=0;j<binsz;j++)
        {
            fprintf(fout, "%12.6lf ", rst+(i+1.0)*delt);
            fprintf(fout, "%12.6lf ", rstz+(j+1.0)*deltz);
            fprintf(fout, "%12.10lf ", his_2D[i][j]);
            fprintf(fout, "\n");           
        }
        fprintf(fout, "\n");   
    }
    fclose(fout);

    sNhB=sNhB*1.0/(natoms /3)/len;
    printf("##### RESULTS #####\n");
    printf("  The len = %d \n", len);
    printf("  N(HB)) = %lf \n", sNhB);
    printf("##### RESULTS #####\n");
    fout=fopen("NBH.dat", "w");
    fprintf(fout, "%lf\n", sNhB);
    fclose(fout);

    printf(" >> Please plot the figures using gnuplot\n");

   
}

