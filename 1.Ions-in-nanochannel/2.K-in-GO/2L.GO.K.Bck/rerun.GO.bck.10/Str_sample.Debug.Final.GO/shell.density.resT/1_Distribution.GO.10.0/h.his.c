#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define N_MOL 96
#define MAX_T 20001
#define MAXL 2000
//const int N_MOL=52;
//const int MAX_T=20001;
// The reason for this warning is that const in c doesn't mean constant. 
// It means "read only". So the value is stored at a memory address and 
// could potentially be changed by machine code.
double MAX(double x, double y)
{
	if (x>y) return x;
	return y;
}
double MIN(double x, double y)
{
	if (x>y) return y;
	return x;
}
void LoadMass(int n, double *mass)
{
	FILE *fin;
	int id, i;
	double m;
	fin=fopen("mass.data", "r");
	for (i=0;i<n;i++)
	{	
		fscanf(fin, "%d%lf", &id, &m);
		mass[7+i]=m;
	}
	fclose(fin);
}


int main()
{
	int  i,j;
	int natoms;
	int N1, N2;
	int id, mol, type;
	double q, rx,ry,rz;
	double mass[20];

	double st, ed, dh, leng, delt; 
	int step, kk;
	double his[MAXL][6]; 
		// his[][0] .. particle O
		// his[][1] .. mass all
		// his[][2] .. particle H
		// his[][3] .. 
		// his[][4] .. 
		// his[][5] .. 
	
	double ss[6];
	int bins;
	char header_char[300];


	FILE *fin;
	FILE *fout;

	LoadMass(2, mass);
	dh=7.0;
	st=-dh/2.0;
	ed=st+dh;
	bins=100;
	delt=(ed-st)/bins;

	for (i=0;i<bins;i++) 
	{
		for (j=0;j<6;j++)
		{
			his[i][j]=0.0;
		}
	}
	step=0;

	printf("Start .....\n");

	fin=fopen("../water.unwrap.xyz", "r");
	
	while (feof(fin)==0)
	{
		step++;
		if (step % 100==0) printf("step = %d\n", step);
		
		for (i=0;i<3;i++) fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &natoms);
		for (i=0;i<5;i++) fgets(header_char, 300, fin);

		for (i=0;i<natoms;i++)
		{
			fscanf(fin, "%d%d", &id, &type);
			fscanf(fin, "%lf%lf%lf\n", &rx, &ry, &rz);
			if (rz>dh/2) rz=rz-dh;
			if (rz<-dh/2) rz=rz+dh;
			if (rz>dh/2 || rz<-dh/2)
			{
				printf("ERROR !!\n");
				return 0;
			}
			
			kk=(rz-st)/delt;
			// for particle O distribution 
			if (type==8) 
			{
				his[kk][0]+=1.0;
			}
			// for particle all distribution 
			his[kk][2]+=1.0;

			// for mass 
			his[kk][1]+=mass[type];
		}
	
	}
	fclose(fin);
	// normalization 
	// the area is zero 
	for (j=0;j<6;j++) ss[j]=0.0;

	for (i=0;i<bins;i++) 
	{
		for (j=0;j<6;j++)
		{
			his[i][j]/=step;
			ss[j]+=his[i][j]*delt;
		}
	}
	for (i=0;i<bins;i++) 
	{
		for (j=0;j<6;j++)
		{
			his[i][j]/=ss[j];
		}
	}

	// sysmetric operation 
	int contactpart;
	double sysave;
	for (i=0;i<bins/2;i++)
	{
		for (j=0;j<6;j++)
		{
			contactpart=bins-1-i;
			sysave=(his[i][j]+his[contactpart][j])/2.0;
			his[i][j]=sysave;
			his[contactpart][j]=sysave;
		}
	}

	printf("St = %12.6lf\n ed = %12.6lf\n ", st, ed);
	printf("\n");
	printf("Total step is %d\n", step);

	fout=fopen("Water.distribution.data", "w");
	fprintf(fout, "## bin = %d\n", bins);

	for (i=0;i<bins;i++)
	{
		fprintf(fout, "%12.6lf ", st+(i+0.6)*delt);
		for (j=0;j<6;j++)
		fprintf(fout, "%12.6lf ", his[i][j]);
		fprintf(fout, "\n");
	}

	fclose(fout);

	return 0;
}
