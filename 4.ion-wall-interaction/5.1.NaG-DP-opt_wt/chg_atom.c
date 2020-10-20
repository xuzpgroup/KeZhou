// Find Dirac point 
// 1. Find the local minimal in [-2, 0]
// 2. Fit the 3 points by y=ax^2+bx+c
// 3. Get the fitting minimal, x_{min}=-b/2a

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define N_C  32
#define N_WT 1
#define MAX_N 100


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



int main()
{
	int len, i, j, ind;
	double r[MAX_N][3], q[MAX_N], dis[MAX_N], vol[MAX_N];
	double s_q;
	char header_char[300];
	char filename[300];

	FILE *fin;
	FILE *fout;


	sprintf(filename, "./ACF.dat");
	fin=fopen(filename, "r");
	printf("Cal .....\n");
	fgets(header_char, 300, fin);
	fgets(header_char, 300, fin);

	// for C atom
	j=0;
	s_q=0.0;
	for (i=0;i<N_C;i++)
	{
		fscanf(fin, "%d", &ind);
		fscanf(fin, "%lf%lf%lf", &r[j][0], &r[j][1], &r[j][2]);
		fscanf(fin, "%lf", &q[j]);
		fscanf(fin, "%lf", &dis[j]);
		fscanf(fin, "%lf", &vol[j]);
		fscanf(fin, "\n");
		s_q+=q[j];
		j++;		
	}
	printf("The Total charge on C is %lf, ave = %lf , \
		rezidue = %lf\n", s_q, s_q/N_C, s_q-4*N_C);

	// for metal atom
	// fgets(header_char, 300, fin);
	s_q=0.0;
		fscanf(fin, "%d", &ind);
		fscanf(fin, "%lf%lf%lf", &r[j][0], &r[j][1], &r[j][2]);
		fscanf(fin, "%lf", &q[j]);
		fscanf(fin, "%lf", &dis[j]);
		fscanf(fin, "%lf", &vol[j]);
		fscanf(fin, "\n");
		s_q+=q[j];
		j++;
	printf("The Total charge on M is %lf, ave = %lf , \
		rezidue = %lf\n", s_q, s_q, s_q-7);

	// for O atom
	s_q=0.0;
	for (i=0;i<N_WT;i++)
	{
		fscanf(fin, "%d", &ind);
		fscanf(fin, "%lf%lf%lf", &r[j][0], &r[j][1], &r[j][2]);
		fscanf(fin, "%lf", &q[j]);
		fscanf(fin, "%lf", &dis[j]);
		fscanf(fin, "%lf", &vol[j]);
		fscanf(fin, "\n");
		s_q+=q[j];
		j++;		
	}
	printf("The Total charge on O is %lf, ave = %lf , \
		rezidue = %lf\n", s_q, s_q/N_WT, s_q/N_WT-6);

	// for H atoms
	s_q=0.0;
	for (i=0;i<N_WT*2;i++)
	{
		fscanf(fin, "%d", &ind);
		fscanf(fin, "%lf%lf%lf", &r[j][0], &r[j][1], &r[j][2]);
		fscanf(fin, "%lf", &q[j]);
		fscanf(fin, "%lf", &dis[j]);
		fscanf(fin, "%lf", &vol[j]);
		fscanf(fin, "\n");
		s_q+=q[j];
		j++;		
	}
	printf("The Total charge on H is %lf, ave = %lf , \
		rezidue = %lf\n", s_q, s_q/N_WT/2, s_q/N_WT/2-1);

	fclose(fin);

	return 0;
}



