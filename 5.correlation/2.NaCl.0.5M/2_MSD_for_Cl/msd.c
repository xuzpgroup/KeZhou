// Calculate MSD curve and fit the diffusion coefficent by Einstan relationship
// MSD calculation for Cl-
// 2017.10.19	the order of atom is random in dumm file
// sort accoeding to atomic ID
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define N_MOL 15
#define MAX_T 100000
#define	DDT 50
//const int N_MOL=52;
//const int MAX_T=20001;
// The reason for this warning is that const in c doesn't mean constant. 
// It means "read only". So the value is stored at a memory address and 
// could potentially be changed by machine code.

double x[MAX_T][N_MOL];
double y[MAX_T][N_MOL];
double z[MAX_T][N_MOL];
//double msd[MAX_T][N_MOL][4];
double xt[N_MOL],yt[N_MOL],zt[N_MOL];


void LinearFitting(double *x, double *y, int n, double stf)
{
	FILE *ffit;
	int i, st, ed;
	double a, b;
	double sx2, sy, sx, sxy;

	ffit=fopen("Fitting.log", "w");
	st=floor(n*stf+0.01);
	ed=n-1;
	sx2=0;
	sy=0;
	sx=0;
	sxy=0;
	printf("The total length of data =%d\n", n);
	printf("Start from  =%d\n", st);
	
	// y=a+bx
	for (i=st;i<=ed;i++)
	{
		sx2+=x[i]*x[i];
		sx+=x[i];
		sy+=y[i];
		sxy+=x[i]*y[i];
	}
	a=(sx2*sy-sx*sxy)/((n-st)*sx2-sx*sx);
	b=((n-st)*sxy-sx*sy)/((n-st)*sx2-sx*sx);

	printf("The fitting results: \n");
	printf("\t\t D=1/6*d(MSD)/dt\n");
	printf("\t\tD = %lf 1.0E-9 m2/s\n", b*10/4.0);
	for (i=st;i<=ed;i++)
	{
		fprintf(ffit, "%lf ", x[i]);
		fprintf(ffit, "%lf ", a+b*x[i]);
		fprintf(ffit, "\n");
	}
	fclose(ffit);
	printf("Finsh fitting\n");
}


int main()
{
	int len, tt; 
	int dt;
	int i,j,st;
	int dem;
	char header_char[300],tp;
	double s[N_MOL][4];
	double s1, s2;
	double timestep;
	int id, mol, type, rid;
	double  q, x1, y1, z1, xu, yu, zu, vx, vy, vz, chg; 
	int kk, limt;
	double fx[1000], fy[1000];

	FILE *dumpxu;
	FILE *fmsd;
	FILE *fdp, *fz;
	FILE *ftar;
	printf("Start .....\n");

	len=0;
	dumpxu=fopen("../ions.dat", "r");
	fmsd=fopen("msd.data", "w");

	
	tt=0;
	timestep=0.250; // ## ps
	
	while (feof(dumpxu)==0)
	{
		len++;
		tt++;
		// printf("%ld\n", len);
		for (i=0;i<9;i++)
		{
			fgets(header_char, 300, dumpxu);
		}
		
		// fprintf(fdp, "%d\n", 3*N_MOL);
		// fprintf(fdp, "step: %d\n", len*1000);
		j=0;
		for (i=0;i<N_MOL*2;i++)
		{
			// ITEM: ATOMS id mol type xu yu zu vx vy vz
			fscanf(dumpxu, "%d%d%d%lf",&id, &mol, &type, &chg);
			fscanf(dumpxu, "%lf%lf%lf", &vx, &vz, &vz);
			fscanf(dumpxu, "%lf%lf%lf", &xu, &yu, &zu);
			fscanf(dumpxu, "\n");				
			if (type==10)
			{
				x[len-1][j]=xu;
				y[len-1][j]=yu;
				z[len-1][j]=zu;
				j++;
			}
		}
		if (tt<000) len--;
	}	
	len--;
	fclose(dumpxu);


	printf("Calaltion .....\n");
	printf("Total length is %d\n", len);

	limt=floor(1000.0/timestep+0.01);
	// printf("The max_T = %d ps\n", limt);

	kk=0;
	for (dt=DDT;dt<=limt;dt=dt+DDT)
	{
		kk++;
		for (dem=0;dem<4;dem++)
		{
			for (j=0;j<N_MOL;j++)
			{
				//msd[dt][j][dem]=0.0;
				s[j][dem]=0.0;	
			}
		} 
		st=0;
		i=0;
		if ( dt%1000==2 )printf("Now is %-12d\n", dt);

		while (st+dt<=len)
		{
			for (j=0;j<N_MOL;j++)
			{
				xt[j]=(x[st][j]-x[st+dt][j])*(x[st][j]-x[st+dt][j]);
				yt[j]=(y[st][j]-y[st+dt][j])*(y[st][j]-y[st+dt][j]);
				zt[j]=(z[st][j]-z[st+dt][j])*(z[st][j]-z[st+dt][j]);
				s[j][1]=s[j][1]+xt[j];
				s[j][2]=s[j][2]+yt[j];
				s[j][3]=s[j][3]+zt[j];
				s[j][0]=s[j][0]+xt[j]+yt[j];
			}
			st=st+1;
			i++;					
		}

		fprintf(fmsd, "%-15.4lf ", dt*timestep);
		// printf("%-12d\n", dt);
		for (dem=0;dem<4;dem++)
		{
			//s1=0.0;
			s2=0.0;
			for (j=0;j<N_MOL;j++)
			{
				//if (j< N_MOL/2) s1=s1+s[j][dem]/i;
				//if (j>=N_MOL/2) s2=s2+s[j][dem]/i;
				s2=s2+s[j][dem]/i;
			}
			//s1=s1/N_MOL*2.0;
			s2=s2/N_MOL;
			fprintf(fmsd, "%-15.4lf ", s2);
			if (dem==0) fprintf(fmsd, "%-15.4lf ", s2/4.0/(timestep*dt)*10); //unit 10-9 m2/s
			if (dem==0)
			{
				fx[kk-1]=dt*timestep;
				fy[kk-1]=s2;
			}

			//printf("%-15.4lf%-15.4lf\n", s1, s2);	
		}
		fprintf(fmsd, "\n");
		// break;
	}
	
	fclose(fmsd);


	// linear fitting the data 
	LinearFitting(fx, fy, kk, 0.2);

	return 0;
}



