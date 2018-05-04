// #####################################################################
// #                                                                   #
// # This program uses Jarzynski relationship to calculate free energy #
// # profiles from work done in independent SMD simulations            #
// #                                                                   #
// #####################################################################

//  Program written by Zhou Ke, zhoukejustin@163.com
//  Last update 14/05/2017


//  Description: This program uses the Jarzynski relationship to convert a collection of
//  work data compiled from multiple steered molecular dynamics simulations into a free
//  energy profile for the process. 
//  Reference: The Journal of Chenmical Physics 119, 3359(2003)

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "math_zhou.h"


#define MAX_LEN    30
#define KB 1.98722E-3
//  KB = Boltzmann's constant in kcal/mol*K (R, or gas constant)


char    floderheader[20];
int     len;
double  vel;
double  temp; //  TEMP = system temperature


void JE_expential()
{
    int i, j;
    int step;
    double fpull[7], work[MAX_LEN], x, y, pmf;
    FILE *fin[20];
    FILE *fout;
    char header_char[300]; 
    char filepath[100], filename[100];
    int flag;

    printf("\n  Implement Jarzynski's method \n");
    printf("  exp(F(x)/kT) = < exp(W(x)/kT) >\n\n");
    
    for (i=0;i<len;i++)
    {
        sprintf(filepath, "%s.%d/%s", floderheader, i+1, "smd.pmf.data");
        puts(filepath);
        fin[i]=fopen(filepath, "r");
        fgets(header_char, 100, fin[i]);
        fgets(header_char, 100, fin[i]);    
    }

    sprintf(filename, "%s.pmf.exp.data", floderheader);
    fout=fopen(filename, "w");
    fprintf(fout, "#######  x (A) #### PMF (Kcal/mol)###########\n");

    flag=0;

    while (!flag)
    {
        for (i=0;i<len;i++)
        {
            sprintf(filepath, "%s.%d/%s", floderheader, i+1, "smd.pmf.data");
            fscanf(fin[i], "%d", &step);
            for (j=0;j<7;j++)
            {
                fscanf(fin[i], "%lf", &fpull[j]);
            }
            fscanf(fin[i], "\n");
            work[i]=fpull[6];
            x=step*1.0E-6*vel;
        }      
        y=array_bolzman_ave(len, work, KB, temp);
        pmf=-log(y)*KB*temp;
        fprintf(fout, "%12.6lf  %12.6lf\n", x, pmf);
        // if (x>12.0) break;
        flag=feof(fin[0]);  // end or not 
    }

    fclose(fout);
    for (i=0;i<len;i++) fclose(fin[i]);
}

void JE_2rd()
{
    int i, j;
    int step;
    double fpull[7], work[MAX_LEN], x, y1, y2, pmf, pmf_biased, pmf_unbiased;
    FILE *fin[20];
    FILE *fout;
    char header_char[300]; 
    char filepath[100], filename[100];
    int flag;

    printf("\n  Implement Jarzynski's 2nd order expansion \n");
    printf("  F = < W > - (<W^2> - <W>^2)/(2*kb*T) \n\n");

    for (i=0;i<len;i++)
    {
        sprintf(filepath, "%s.%d/%s", floderheader, i+1, "smd.pmf.data");
        puts(filepath);
        fin[i]=fopen(filepath, "r");
        fgets(header_char, 100, fin[i]);
        fgets(header_char, 100, fin[i]);    
    }

    sprintf(filename, "%s.pmf.2nd.data", floderheader);
    fout=fopen(filename, "w");
    fprintf(fout, "#######  x (A) #### PMF (Kcal/mol) biased and unbiased estimation ###########\n");

    flag=0;

    while (!flag)
    {
        for (i=0;i<len;i++)
        {
            sprintf(filepath, "%s.%d/%s", floderheader, i+1, "smd.pmf.data");
            fscanf(fin[i], "%d", &step);
            for (j=0;j<7;j++)
            {
                fscanf(fin[i], "%lf", &fpull[j]);
            }
            fscanf(fin[i], "\n");
            work[i]=fpull[6];
            x=step*1.0E-6*vel;
        }      
        y1=array_ave(len, work);
        y2=array_square_ave(len, work);
        pmf_biased  = y1-(y2-y1*y1)/2/KB/temp;
        pmf_unbiased= y1-(y2-y1*y1)/2/KB/temp*len/(len-1);
        // printf("%d  %d  %lf\n", len, len-1, (len+0.0)/(len-1.0));    
        // printf("<W>= %lf, <W^2>= %lf, <W>^2= %lf  \n", y1, y2, y1*y1);    
        fprintf(fout, "%12.6lf  %12.6lf  %12.6lf\n", x, pmf_biased, pmf_unbiased);
        // if (x>12.0) break;
        flag=feof(fin[0]);  // end or not 
    }

    fclose(fout);
    for (i=0;i<len;i++) fclose(fin[i]);
}

int main()
{
    FILE *f0;
    f0=fopen("PMF.Path.info", "r");
    fscanf(f0, "%s\n", floderheader);
    fscanf(f0, "%d\n", &len);
    fscanf(f0, "%lf\n", &vel);
    fscanf(f0, "%lf\n", &temp);
    fclose(f0);

    printf("#####################################################################\n");
    printf("#                                                                   #\n");
    printf("# This program uses Jarzynski relationship to calculate free energy #\n");
    printf("#                                                                   #\n");
    printf("#####################################################################\n");
    
    
    printf("\n");
    printf("%d Samples \n", len);
    printf("vel = %.2f A/ns \n", vel);
    printf("temp = %.2f K\n", temp);
    printf("\n");
    
    JE_expential();
    JE_2rd();

   
	return 0;
}

