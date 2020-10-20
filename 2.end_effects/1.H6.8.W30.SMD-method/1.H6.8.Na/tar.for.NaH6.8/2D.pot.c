// average basing on histogarm

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define MAXNATOM 20000
#define Rcut 3.3 
#define HisItems 1
#define MAXBINS 500
#define Height 10.0


// for histogram
int    bins[HisItems][2];
double hst[HisItems][2], hed[HisItems][2], delt[HisItems][2];
double his[MAXBINS][MAXBINS][HisItems];
	// 	info
	// 	0	num for O atom
	// 	1	mass
double box[3];

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
double pbc(double ddx, double lx)
{
	double dx;
	dx=ddx;
	while (dx>lx)  dx=dx-lx;
	while (dx<0.0) dx=dx+lx;
	return dx;
}
void PutIntoBox(double *r, double *br, double *box)
{
	br[0]=pbc(r[0], box[0]);
	br[1]=pbc(r[1], box[1]);
	br[2]=pbc(r[2], box[2]);
}


void HistogramIt(double samples1, double samples2, int jj, double factor)
{
	int kk1, kk2; 
	kk1=(samples1-hst[jj][0])/delt[jj][0];
	kk2=(samples2-hst[jj][1])/delt[jj][1];
	// printf("%d %d \n", kk1, kk2);
	if (equal(samples1, hed[jj][0])) kk1--;
	if (equal(samples2, hed[jj][1])) kk2--;
	
	
	// printf("%lf %d  jj=%d \n", samples, kk, jj);
	if (kk1<bins[jj][0] && kk1>=0 && kk2<bins[jj][1] && kk2>=0)
	{
		his[kk1][kk2][jj]+=1.0*factor;
		// other 
	}
	else {
		printf("---%lf %lf %d %d \n", samples1, samples2, kk1, kk2);
		printf("Error when histogram! \n");
	}
}

int main()
{
	int len;
	int step, natoms, nwater;
	int i, j, l, st, tp, id;
	double r0[3], ri[3], ti, pw;

	char header_char[300],tp1,tp2;
	char filename[300];
	double data[HisItems];
		

	FILE *fin;
	FILE *fout;

	box[0]=8.4000000000000000e+01;
	box[1]=3.3948000000000000e+01;
	box[2]=6.8*2;
	pw=52.4924;
	
	sprintf(filename, "./pos.data");
	
	
//	histogram and ave
	printf("Start .....\n");

	for (i=0;i<HisItems;i++) bins[i][0]=100;
	for (i=0;i<HisItems;i++) bins[i][1]=100;

	for (i=0;i<HisItems;i++)
	{
		hst[i][0]=40;		// min  ... x direction 
		hed[i][0]=64;	// max
		hst[i][1]=0;		// min  ... z direction 
		hed[i][1]=box[2];	// max  		
		delt[i][0]=(hed[i][0]-hst[i][0])/bins[i][0];
		delt[i][1]=(hed[i][1]-hst[i][1])/bins[i][1];
			
	}
	printf("Delt x =%lf\n", delt[0][0]);
	printf("Delt z =%lf\n", delt[0][1]);


	for (i=0;i<HisItems;i++)
	{
		for (j=0;j<bins[i][0];j++) //... x direction 
		{
			for (l=0;l<bins[i][1];l++) // ... y direction 
			{
				his[j][l][i]=0;
				// other 
			}
		}
	}

	fin=fopen(filename, "r");
	printf("Cal .....\n");
	// fgets(header_char, 300, fin);
	// fgets(header_char, 300, fin);

	len=0;
	while (feof(fin)==0)
	{
		len++;
		if (len % 1000 ==0) printf("Now is %d\n", len);
		// input raw data
		fscanf(fin, "%lf%lf%lf%lf\n", &ti, &r0[0], &r0[1] ,&r0[2]);
		HistogramIt(r0[0], r0[2], 0, 1.0);
	}	
	fclose(fin);

	// time average the data 
	// for (i=0;i<HisItems;i++)
	// {
	// 	for (j=0;j<bins[i][0];j++)
	// 	{
	// 		for (l=0;l<bins[i][1];l++)
	// 		{
	// 			his[j][l][i]/=len;
	// 			// other 
	// 		}
	// 	}
	// }


	//output 
	fout=fopen("./His.2D.ions.dat", "w");
	fprintf(fout, "## bin x bin = %d x %d\n", bins[1][0], bins[1][1]);
	printf("output ..\n Total len = %d\n", len);
	int PBC, ipbc;
	PBC=bins[0][0]/bins[0][1];


	for (i=0;i<bins[0][1];i++)
	{
		for (ipbc=0;ipbc<PBC;ipbc++)
		{
			for (j=0;j<bins[0][0];j++)
			{
				for (l=0;l<HisItems;l++)
				{
					// printf("%d %d %d %d \n", i, j, ipbc, l);
					fprintf(fout, "%12.6lf ", hst[l][0]+(j+0.5)*delt[l][0]);
					fprintf(fout, "%12.6lf ", hst[l][1]+(i+0.5)*delt[l][1]+ipbc*box[2]);
					fprintf(fout, "%12.6f ",his[j][i][l]);
				}	
				fprintf(fout, "\n");
			}
		}
		fprintf(fout, "\n");
	}
	fclose(fout);


	//output  1D 
	printf("output 1D info \n");
	double ss;
	fout=fopen("./His.1D.dat", "w");
	fprintf(fout, "## bin = %d \n", bins[0][1]);
	printf("output ..\n Total len = %d\n", len);
	for (i=0;i<bins[0][0];i++)
	{
		fprintf(fout, "%12.6lf ", hst[0][0]+(i+0.5)*delt[0][0]);
		ss=0.0;
		for (j=0;j<=bins[0][1];j++)
		{
			ss+=his[i][j][0];	
		}
		fprintf(fout, "%12.6lf ", ss);	
		fprintf(fout, "\n");
	}

	fclose(fout);


	return 0;
}



