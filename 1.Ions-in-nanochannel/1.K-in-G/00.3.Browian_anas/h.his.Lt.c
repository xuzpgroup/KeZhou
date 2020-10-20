#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define MAXNATOM 20000
#define Rcut 3.3 
#define HisItems 3
#define MAXBINS 500


// for histogram
int    bins;
double hst[HisItems], hed[HisItems], delt[HisItems];
double his[MAXBINS][HisItems];
	// 	info
	// 	0	Lt_x
	// 	1	Lt_y
	// 	2 	Lt_z

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


void HistogramIt(double samples, int jj, double factor)
{
	int kk; 
	kk=(samples-hst[jj])/delt[jj];
	// printf("%lf %d  jj=%d \n", samples, kk, jj);
	if (kk<bins)
	{
		his[kk][jj]+=1.0*factor;
	}
}

int main()
{
	int len;
	int step, natoms, nwater;
	int i, j, st, tp;

	char header_char[300];
	char filename[300];
	double data[HisItems];
	double range[HisItems][2];	// max and min
	double ss[HisItems];
		

	FILE *fin;
	FILE *fout;


	// get the range of raw data 
	for (i=0;i<HisItems;i++)
	{
		range[i][0]=0.0; // min
		range[i][1]=0.0; // max	
	}
	fin=fopen("Lt.tot.Lt.20.dat", "r");
	while (feof(fin)==0)
	{
		for (i=0;i<HisItems;i++) fscanf(fin, "%lf%lf%lf\n", &data[0], &data[1], &data[2]);
		for (i=0;i<HisItems;i++)
		{
			range[i][1]=DMax(range[i][1], DABS(data[i]));
		}	
	}
	fclose(fin);
	for (i=0;i<HisItems;i++)
	{
		printf("Min %d = %lf; Max %d = %lf\n", i+1, range[i][0], i+1, range[i][1]);
	}

//	histogram
	printf("Start .....\n");

	bins=50;
	for (i=0;i<HisItems;i++)
	{
		hst[i]=range[i][0];	// min
		hed[i]=range[i][1];	// max
		delt[i]=(hed[i]-hst[i])/bins;
	}

	for (i=0;i<bins;i++)
	{
		for (j=0;j<HisItems;j++)
			his[i][j]=0.0;
	}

	fin=fopen("Lt.tot.Lt.20.dat", "r");
	printf("Cal .....\n");

	len=0;
	while (feof(fin)==0)
	{
		len++;
		if (len % 10000 ==0) printf("Now is %d\n", len);
		// input raw data
		for (i=0;i<HisItems;i++) fscanf(fin, "%lf%lf%lf\n", &data[0], &data[1], &data[2]);
		HistogramIt(DABS(data[0]), 0, 1.0);		// Ltx
		HistogramIt(DABS(data[1]), 1, 1.0);		// Lty
		HistogramIt(DABS(data[2]), 2, 1.0);		// Ltz
	}	
	fclose(fin);

	// normalization 
	// the area is zero 
	for (j=0;j<HisItems;j++) ss[j]=0.0;
	for (i=0;i<bins;i++) 
	{
		for (j=0;j<HisItems;j++)
		{
			// his[i][j]/=len;
			ss[j]+=his[i][j]*delt[j];
		}
	}
	for (i=0;i<bins;i++) 
	{
		for (j=0;j<HisItems;j++)
		{
			his[i][j]/=ss[j];
		}
	}	
	//output 
	fout=fopen("./His.Lt.20ps.dat", "w");
	fprintf(fout, "## bin = %d\n", bins);
	for (i=0;i<bins;i++)
	{
		for (j=0;j<HisItems;j++)
		{
			fprintf(fout, "%12.6lf ", hst[j]+(i+0.5)*delt[j]);
			fprintf(fout, "%12.6lf ", his[i][j]);
		}
		fprintf(fout, "\n");
	}
	fclose(fout);

	return 0;
}



