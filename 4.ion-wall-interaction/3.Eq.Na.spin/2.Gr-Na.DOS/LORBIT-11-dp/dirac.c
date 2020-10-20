// Find Dirac point 
// 1. Find the local minimal in [-2, 0]
// 2. Fit the 3 points by y=ax^2+bx+c
// 3. Get the fitting minimal, x_{min}=-b/2a

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define MAXLEN 500


double DABS(double x)
{
	if (x<0) return -x;
	return x;
}
double DMax(double x, double y)
{
	if (x>y) return x;
	return y;
}

int equal(double x, double y)
{
	if (DABS(x-y)<0.01) return 1;
	return 0;
}

double Fitting_3p_get_minimal(double *x, double *y)
{
	double s1, s2;
	double dx21, dx22;
	double dx1, dx2;
	double dy1, dy2;

	dx1=x[0]-x[1];
	dx2=x[1]-x[2];
	dy1=y[0]-y[1];
	dy2=y[1]-y[2];	
	dx21=x[0]*x[0]-x[1]*x[1];
	dx22=x[1]*x[1]-x[2]*x[2];

	s2=dy1*dx2-dy2*dx1;
	s1=dx21*dy2-dx22*dy1;

	return -s1/2.0/s2;
}

int main()
{
	int len, j;
	double mine, min1, min2;
	double ener[MAXLEN], sp_up[MAXLEN], sp_dw[MAXLEN];
	double x[3], y[3];
	char header_char[300];
	char filename[300];

		

	FILE *fin;
	FILE *fout;


	sprintf(filename, "./TDOS.dat");
	fin=fopen(filename, "r");
	printf("Cal .....\n");


	len=0;
	mine=100.0;
	j=0;

	fgets(header_char, 300, fin);
	while (feof(fin)==0)
	{
		// input raw data
		len++;
		fscanf(fin, "%lf", &ener[len-1]);
		fscanf(fin, "%lf%lf", &sp_up[len-1], &sp_dw[len-1]);
		sp_dw[len-1]=sp_dw[len-1]*-1;
		fscanf(fin, "\n");

		if (ener[len-1]<0.0 && ener[len-1]>-2.0)
		{
			if (sp_up[len-1]+sp_dw[len-1]<mine*2)
			{
				j=len-1;
				mine=0.5*(sp_up[len-1]+sp_dw[len-1]);
			}			
		}
		
	}	
	printf("Total length = %d\n", len);
	printf("The min point is %lf with min_e %lf (%lf, %lf)\n", ener[j], mine,\
		sp_up[j], sp_dw[j]);
	fclose(fin);
//  fitting 
	printf("Get Dirac point by fitting DOS .. \n");
	x[1]=ener[j];
	x[0]=ener[j-1];
	x[2]=ener[j+1];

	y[1]=sp_up[j];
	y[0]=sp_up[j-1];
	y[2]=sp_up[j+1];
	min1=Fitting_3p_get_minimal(x, y);

	y[1]=sp_dw[j];
	y[0]=sp_dw[j-1];
	y[2]=sp_dw[j+1];
	min2=Fitting_3p_get_minimal(x, y);
	
	printf("The min point is min_1 = %lf  min_2 =%lf,\n average = %lf\n",\
		min1, min2, (min1+min2)/2.0);	

	sprintf(filename, "./DIRAC.dat");
	fout=fopen(filename, "w");
	fprintf(fout, "%lf\n", (min1+min2)/2.0);
	fprintf(fout, "%lf %lf \n", min1, min2);
	fclose(fout);

	return 0;
}



