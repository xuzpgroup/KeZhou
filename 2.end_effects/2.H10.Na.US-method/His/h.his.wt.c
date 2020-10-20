#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define N_MOL 96
#define MAX_T 20001
#define MAXL 2000
#define NSIS 1
#define SAMPLES 26
#define	M_middle 2.55

double Zsample[SAMPLES];

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
double DABS(double x)
{
	if (x<0) return -x;
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

void ALL_sample()
{
	// Zsample[0]=1.2;	
	// Zsample[1]=1.1;
	// Zsample[2]=1.0;
	// Zsample[3]=0.9;
	// Zsample[4]=0.8;
	// Zsample[5]=0.7;
	// Zsample[6]=0.6;
	// Zsample[7]=0.5;
	// Zsample[8]=0.4;
	// Zsample[9]=0.3;
	// Zsample[10]=0.2;
	// Zsample[11]=0.1;
	// Zsample[12]=0.0;
	// Zsample[13]=-0.1;
	// Zsample[14]=-0.2;
	// Zsample[15]=-0.3;
	// Zsample[16]=-0.4;
	// Zsample[17]=-0.5;
	// Zsample[18]=-0.6;
	// Zsample[19]=-0.7;
	// Zsample[20]=-0.8;
	// Zsample[21]=-0.9;
	// Zsample[22]=-1.0;
	// Zsample[23]=1.3;
	// Zsample[24]=1.4;
	// Zsample[25]=1.5;
	// Zsample[26]=1.6;
	// Zsample[27]=1.7;
	// Zsample[28]=1.8;
	// Zsample[29]=1.9;
	// Zsample[30]=2.0;
	// Zsample[31]=2.1;
	// Zsample[32]=2.2;
	// Zsample[33]=2.3;
	// Zsample[34]=2.4;
	// Zsample[35]=2.5;
	Zsample[0]=1.0;
	Zsample[1]=0.9;
	Zsample[2]=0.8;
	Zsample[3]=0.7;
	Zsample[4]=0.6;
	Zsample[5]=0.5;
	Zsample[6]=0.4;
	Zsample[7]=0.3;
	Zsample[8]=0.2;
	Zsample[9]=0.1;
	Zsample[10]=0.0;
	Zsample[11]=-0.1;
	Zsample[12]=-0.2;
	Zsample[13]=-0.3;
	Zsample[14]=-0.4;
	Zsample[15]=-0.5;
	Zsample[16]=-0.6;
	Zsample[17]=-0.7;
	Zsample[18]=-0.8;
	Zsample[19]=-0.9;
	Zsample[20]=-1.0;		
	Zsample[21]=1.1;
	Zsample[22]=1.2;
	Zsample[23]=1.3;
	Zsample[24]=1.4;
	Zsample[25]=1.5;

}


int main()
{
	int  i,j, ii;
	int natoms;
	int N1, N2;
	int id, mol, type;
	double  q, x1, y1, z1, xu, yu, zu; 
	double rx,ry,rz;
	double mass[20];
	double da[10];
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
	char filename[300];

	FILE *fin;
	FILE *fout;
	FILE *fsample;

	// LoadMass(2, mass);
	ALL_sample();

	st=-2.0;
	ed=3.0;
	bins=500;
	delt=(ed-st)/bins;


	step=0;

	fout=fopen("distribution.data", "w");

	printf("Start .....\n");
	for (ii=0;ii<SAMPLES;ii=ii+1)
	{
		for (i=0;i<bins;i++) 
		{
			for (j=0;j<NSIS;j++)
			{
				his[i][j]=0.0;
			}
		}
		// if (ii==6) continue;
		// if (ii>=12) sprintf(filename, "../COLVAR.%.2lf", Zsample[ii]);
		// if (ii<=11) sprintf(filename, "../COLVAR.10ns.%.2lf", Zsample[ii]);
		sprintf(filename, "../COLVAR.10ns.%.2lf", Zsample[ii]);
		puts(filename);
		fin=fopen(filename, "r");
		// sprintf(filename, "./COLVAR.rev.%.2lf", M_middle-Zsample[ii]);
		sprintf(filename, "./COLVAR.rev.%.2lf", Zsample[ii]);
		fsample=fopen(filename, "w");
		fgets(header_char, 300, fin);
		while (feof(fin)==0)
		{
			fscanf(fin, "%lf%lf%lf\n", &da[0], &da[1], &da[2]);
			// if (da[0]>9000) continue;
			// da[1]=DABS(da[1]);
			da[1]=1.0*da[1];
			fprintf(fsample, "%lf  %lf\n", da[0], da[1]);
			if (da[1]<st || da[1]>ed) continue;
			kk=(da[1]-st)/delt;
			his[kk][0]+=1.0;
		}
		fclose(fin);
		fclose(fsample);
		// normalization 

		sprintf(filename, "#### COLVAR.%.1lf\n", Zsample[ii]);
		// sprintf(filename, "#### COLVAR.%.2lf\n", M_middle-Zsample[ii]);
		fputs(filename, fout);
		for (i=0;i<bins;i++)
		{
			fprintf(fout, "%12.6lf ", st+(i+0.5)*delt);
			for (j=0;j<NSIS;j++)
			fprintf(fout, "%12.6lf ", his[i][j]);
			fprintf(fout, "\n");
		}
		fprintf(fout, "\n");

	}



	return 0;
}
