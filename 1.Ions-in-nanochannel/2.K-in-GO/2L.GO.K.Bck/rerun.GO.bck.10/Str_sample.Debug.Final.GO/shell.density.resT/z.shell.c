#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define MAXNATOM 20000
#define Rcut 3.3 
#define HisItems 2
#define MAXBINS 500
#define Height 10.2


// for histogram
int    bins;
double hst[HisItems], hed[HisItems], delt[HisItems];
double his[MAXBINS][HisItems];
	// 	info
	// 	0	O atoms
	// 	1   mass

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

int SamePoint(double *x, double *y)
{
	if (!equal(x[0], y[0])) return 0;
	if (!equal(x[1], y[1])) return 0;
	if (!equal(x[2], y[2])) return 0;
	return 1;
}


void HistogramIt(double samples, int jj, double factor)
{
	int kk; 
	kk=(samples-hst[jj])/delt[jj];
	// printf("%lf %d  jj=%d \n", samples, kk, jj);
	if (kk<bins && kk>=0)
	{
		his[kk][jj]+=1.0*factor;
	}
	else {
		printf("Error when histogram! \n");
	}
}

int main()
{
	int len;
	int natoms, tp;
	int i, j;
	double x, y, z;
	char header_char[300];
	char filename[300];
	double data[HisItems];
	double range[HisItems][2];	// max and min
	double ss[HisItems];
	double mass[20];		
	double timestep;
	
	// for staying time 
	int shell_list[10];
	double shell_r[10][3];
	
	int list_now[10];
	double shell_now[10][3];
	
	int flagOH, lastflagOH;
	int OHtime, OHsep;


	FILE *fin;
	FILE *fout;

	sprintf(filename, "../shell.str.xyz");
	// get the range of raw data 
	for (i=0;i<HisItems;i++)
	{
		range[i][0]=-Height/2.0; // min
		range[i][1]=Height/2.0; // max	
	}

	mass[7]=1.0;	// H atom
	mass[8]=16.0;	// O atom
	timestep=0.05;
	
	
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

	fin=fopen(filename, "r");
	printf("Cal .....\n");

	shell_list[0]=0;
	flagOH=0;
	lastflagOH=0;
	list_now[0]=0;
	
	OHtime=0;
	OHsep=0;


	len=0;
	while (feof(fin)==0)
	{
		len++;
		flagOH=0;
		list_now[0]=0;
		// if (len==10) break;
		if (len %1000==0) printf("Now is %d\n", len);
		fscanf(fin, "%d\n", &natoms);
		fgets(header_char, 300, fin);
		fgets(header_char, 300, fin); // for ions 
		for (i=1;i<natoms;i++)
		{
			fscanf(fin, "%d", &tp);
			fscanf(fin, "%lf%lf%lf", &x, &y, &z);
			// printf("%d %lf %lf %lf \n", tp, x, y, z);
			fscanf(fin, "\n");
			if (tp==8) //O atom 
			{
				HistogramIt(z, 0, 1.0);
				HistogramIt(z, 1, mass[8]);			
			}
			if (tp==7) // H atom
			{
				HistogramIt(z, 1, mass[7]);			
			}
			if (tp==4)
			{
				list_now[0]++;
				shell_now[list_now[0]][0]=x;
				shell_now[list_now[0]][1]=y;
				shell_now[list_now[0]][2]=z;
				flagOH=1;					
			}
		}
		if (flagOH) OHtime++;
		if (flagOH && !lastflagOH)
		{
			OHsep++;
			// printf("gen %d |", len);
		}
		if (!flagOH && lastflagOH)
		{
			// printf("|die %d\n", len);
		}
		// printf("%d %d\n", OHtime, OHsep);
		lastflagOH=flagOH;

		// printf("------------\n");
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

	printf("%d\n", OHtime);
	printf("%d\n", OHsep);
	printf("step= %lf  t= %lf\n", (OHtime+0.0)/OHsep, (OHtime+0.0)/OHsep*timestep);
	fout=fopen("./His.shell.z.dat", "w");
	fprintf(fout, "## bin = %d\n", bins);

	double sysave;

	for (i=0;i<bins;i++)
	{
		for (j=0;j<HisItems;j++)
		{
			sysave=(his[i][j]+his[bins-1-i][j])/2.0;
			fprintf(fout, "%12.6lf ", hst[j]+(i+0.5)*delt[j]);
			// fprintf(fout, "%12.6lf ", his[i][j]);
			fprintf(fout, "%lf", sysave);
		}
		fprintf(fout, "\n");
	}
	fclose(fout);

	return 0;
}



