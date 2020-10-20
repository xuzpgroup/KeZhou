#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_ATOMS 400
#define N_PLANE   60
#define N_WATER	  3
#define  _CC      1.4000

#define SQRT3	  1.732050808

#define BINS 32

double His[BINS][BINS];
double  st[2];
double	ed[2];
double	delt[2];

double MAX_D(double x, double y)
{
    if (x>y) return x;
    return y;
}
double MIN_D(double x, double y)
{
    if (x<y) return x;
    return y;
}

void UNCELL(double *r, double *orgin, double *ur)
{
	double x, y;
	double nx, ny;
	x=r[0]-orgin[0];
	y=r[1]-orgin[1];	
	
	// ny=2.0*y/3.0/_CC;
	// nx=(x-SQRT3/3*y)/(SQRT3*_CC);
	/*for cubic cell*/
	// nx=x/3.0/_CC;
	// ny=y/SQRT3/_CC;
	/*for XX cell*/
	nx=0.5*(x/1.50/_CC-y/(SQRT3/2.0)/_CC);
	ny=0.5*(x/1.50/_CC+y/(SQRT3/2.0)/_CC);
	
	
	nx=nx-floor(nx);
	ny=ny-floor(ny);

	ur[0]=nx;
	ur[1]=ny;
	ur[2]=0.0;
}

void U2XYZ(double *ur, double *rr)
{
	// rr[0]=SQRT3*_CC*ur[0]+SQRT3/2*_CC*ur[1];
	// rr[1]=3.0/2.0*_CC*ur[1];
	/*for cubic cell*/
	// rr[0]=3.0*_CC*ur[0];
	// rr[1]=SQRT3*_CC*ur[1];
	// rr[2]=0.0;
	/*for XX cell*/
	rr[0]=1.5*_CC*(ur[0]+ur[1]);
	rr[1]=SQRT3/2.0*_CC*(-ur[0]+ur[1]);
	
}



void HistogramIt(double samples1, double samples2, double factor)
{
	int kk1, kk2; 
	kk1=(samples1-st[0])/delt[0];
	kk2=(samples2-st[1])/delt[1];
	
	// printf("%lf %d  jj=%d \n", samples, kk, jj);
	if (kk1<=BINS && kk1>=0 && kk2<BINS && kk2>=0)
	{
		His[kk1][kk2]+=1.0*factor;
		// other 
	}
	else {
		printf("---%lf %lf %d %d \n", samples1, samples2, kk1, kk2);
		printf("Error when histogram! \n");
	}
}


int main()
{
	int 	i, j, len;
	char 	header_char[300];
	char 	filename[300], tp;
	int		natoms;
	double 	data[20];
	double  RIons[3];
	double  Orgin[3], nr[3], rr[3];
	double  s1, s2;
	FILE 	*fin, *fout;
	double MINF, MAXF;	




	st[0]=0;
	st[1]=0;
	ed[0]=1.0;
	ed[1]=1.0;
	for (i=0;i<2;i++) delt[i]=(ed[i]-st[i])/BINS;

	Orgin[0]=0.0;
	Orgin[1]=0.00000;
	Orgin[2]=0.00000;

	Orgin[0]=2.09992;
	Orgin[1]=2.4247;
	Orgin[2]=0.00000;

	Orgin[0]=1.3998465;
	Orgin[1]=3.63717;
	Orgin[2]=0.00000;
	
	

	// Orgin[0]=7.274535;
	// Orgin[1]=3.500000;
	// Orgin[2]=7.00000;

	for (i=0;i<BINS;i++)
	{
		for (j=0;j<BINS;j++)
			His[i][j]=0.0;
	}
		
		


	
	fin=fopen("../ions.dat", "r");
	fout=fopen("Pos.dat", "w");

	len=0;
	while (!feof(fin))
	{
		len++;
		if (len % 100000 ==0) printf("Now len = %d\n", len);
		// fscanf(fin, "%d\n", &natoms);
		for (i=0;i<9;i++)
		{
			fgets(header_char, 300, fin);
		}
		for (i=0;i<9;i++) fscanf(fin, "%lf", &data[i]);
		fscanf(fin, "\n");
		RIons[0]=data[3];
		RIons[1]=data[4];
		RIons[2]=data[5];

		UNCELL(RIons, Orgin, nr);
		U2XYZ(nr, rr);
		fprintf(fout, " %lf %lf %lf ", len*0.5, RIons[0], RIons[1]);
		fprintf(fout, " %lf %lf %lf %lf \n", rr[0], rr[1], nr[0], nr[1]);
		// fprintf(fout, " %lf %lf \n", nr[0], nr[1]);
		// HistogramIt(rr[0], rr[1], 1.0);
		HistogramIt(nr[0], nr[1], 1.0);
	}
	printf("TOTAL length = %d \n", len);
	printf("TOTAL Time   = %.2lf ns\n", len*10E-6);

	fclose(fin);
	fclose(fout);

	MINF=1000.0;
	MAXF=-1000;
	for (i=0;i<BINS;i++)
	{
		for (j=0;j<BINS;j++)
		{
			if (His[i][j]<1.0) 
			{
				printf("ERROR %d %d %lf\n", i, j, His[i][j]);
				His[i][j]=1.0;
			}
			
			MINF=MIN_D(His[i][j], MINF);
			MAXF=MAX_D(His[i][j], MAXF);	
		}
	}
	printf(" Min = %lf, Max = %lf\n delt F = %lf kBT\n", log(MINF), log(MAXF), log(MAXF)-log(MINF));

	// sysmetry operation 
	double avehis;
	// for (i=0;i<BINS;i++)
	// {
	// 	for (j=i+1;j<BINS;j++)
	// 	{
	// 		avehis=(His[i][j]+His[j][i])/2.0;
	// 		His[i][j]=avehis;
	// 		His[j][i]=avehis;	
	// 	}
	// }
	// for (i=0;i<BINS;i++)
	// {
	// 	for (j=i;j<BINS;j++)
	// 	{
	// 		avehis=(His[i][j]+His[BINS-1-j][BINS-1-i])/2.0;
	// 		His[i][j]=avehis;
	// 		His[BINS-1-j][BINS-1-i]=avehis;	
	// 	}
	// }

	//output 
	fout=fopen("./His.2D.dat", "w");
	fprintf(fout, "## bin x bin = %d x %d\n", BINS, BINS);
	printf("output ..\n Total len = %d\n", len);
	
	int PBC, ix, jy;
	double fr[3], pr[3];
	// PBC=BINS[0]/BINS[1];
	PBC=3;
	for (ix=-PBC;ix<PBC;ix++)
	{
		for (i=0;i<BINS;i++)
		{
				for (jy=-PBC;jy<PBC;jy++)
				{
					for (j=0;j<BINS;j++)
					{
						fr[0]=jy+ st[0]+(j+0.5)*delt[0];
						fr[1]=ix+ st[1]+(i+0.5)*delt[1];
						fr[2]=0.0;
						U2XYZ(fr, pr);
						// fprintf(fout, "%12.6lf ", st[0]+(j+0.5)*delt[0]);
						// fprintf(fout, "%12.6lf ", st[1]+(i+0.5)*delt[1]);
						fprintf(fout, "%12.6lf %12.6lf ", pr[0], pr[1]);
						// fprintf(fout, "%12.6lf %12.6lf ", fr[0], fr[1]);
						if (His[j][i]<1.0) His[j][i]=1.0;
						fprintf(fout, "%12.6f %12.6f", His[j][i], -log(His[j][i])+log(MAXF));
						fprintf(fout, "\n");						
					}

				}
			fprintf(fout, "\n");
		}		
	}

	fclose(fout);
	
	//output minimal path 
	double Maxmin;
	Maxmin=0.0;
	fout=fopen("./His.path.min.dat", "w");
	for (i=0;i<BINS;i++)
	{
		fprintf(fout, "%lf %lf %lf %lf %lf %lf\n", \
			(i+0.5)/BINS, \
			(i+0.5)/BINS*SQRT3*_CC, \
			-log(His[i][0])+log(MAXF), \
			-log(His[0][i])+log(MAXF), \
			-log(His[i][BINS-1-i])+log(MAXF),\
			-(log(His[i][0])+log(His[0][i])+log(His[i][BINS-1-i]))/3.0+log(MAXF));
		Maxmin=MAX_D(Maxmin, -(log(His[i][0])+log(His[0][i])+log(His[i][BINS-1-i]))/3.0+log(MAXF));
	}
	fclose(fout);
	printf("output ..\n Min F in min path = %lf\n kBT", Maxmin);

	//output maxmial path 
	fout=fopen("./His.path.max.dat", "w");
	for (i=0;i<BINS;i++)
	{
		fprintf(fout, "%lf %lf %lf \n", \
			(i+0.5)/BINS, \
			(i+0.5)/BINS*3*_CC, \
			-log(His[i][i])+log(MAXF));
	}
	fclose(fout);


	return 0;
}



