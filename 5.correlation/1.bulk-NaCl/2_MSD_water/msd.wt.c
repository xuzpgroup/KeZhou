#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define MAXN 2000
#define MAX_T 50000
//const int N_MOL=52;
//const int MAX_T=20001;
// The reason for this warning is that const in c doesn't mean constant. 
// It means "read only". So the value is stored at a memory address and 
// could potentially be changed by machine code.

double x[MAX_T][MAXN];
double y[MAX_T][MAXN];
double z[MAX_T][MAXN];
//double msd[MAX_T][N_MOL][4];
double xt[MAXN],yt[MAXN],zt[MAXN];


int main()
{
	int len, tt; 
	int dt;
	int i,j,st;
	int dem ,natoms, n, tp, id;
	char header_char[300];
	double s[MAXN][4];
	double s1, s2, xu, yu, zu;
	double timestep;


	FILE *dumpxu;
	FILE *fmsd;
	printf("Start .....\n");

	len=0;
	dumpxu=fopen("../sys.xyz", "r");
	//dumpxu=fopen("mol.position.xyz", "r");
	fmsd=fopen("msd.data", "w");
	//printf("%ld\n", len);
	//fgets(header_char, 300, dumpxu);
	//fgets(header_char, 300, dumpxu);
	//fgets(header_char, 300, dumpxu);
	
	tt=0;
	timestep=1000*1.0E-3;  // ps 
	printf("loading ...");	
	while (feof(dumpxu)==0)
	{
		// len++;
		tt++;
		// printf("%ld\n", len);
		for (i=0;i<3;i++) fgets(header_char, 300, dumpxu);
		fscanf(dumpxu, "%d\n", &natoms);
		for (i=0;i<5;i++) fgets(header_char, 300, dumpxu);
		n=0;
		j=0;
		for (i=0;i<natoms;i++)
		{
			fscanf(dumpxu, "%d%d%lf%lf%lf\n", &id, &tp, &xu, &yu, &zu);
			if (tp==8)
			{
				x[len][j]=xu;
				y[len][j]=yu;
				z[len][j]=zu;
				n++;		
				j++;		
			}
		}
		// if (tt % 20 == 1) len++;
		len++;
		if (len % 100==0) printf("Loading %d \n", len);
	}	
	fclose(dumpxu);

	printf("Total length is %d\n", tt);
	printf("Selected length is %d\n", len);
	natoms=n;


	printf("Calaltion .....\n");
	for (dt=1;dt<=1000;dt=dt+1)
	{
		for (dem=0;dem<4;dem++)
		{
			for (j=0;j<natoms;j++)
			{
				//msd[dt][j][dem]=0.0;
				s[j][dem]=0.0;	
			}
		} 
		st=0;
		i=0;
		if ( dt%100==1 )printf("Now is %-12d\n", dt);

		while (st+dt<=len-1)
		{
			for (j=0;j<natoms;j++)
			{
				xt[j]=(x[st][j]-x[st+dt][j])*(x[st][j]-x[st+dt][j]);
				yt[j]=(y[st][j]-y[st+dt][j])*(y[st][j]-y[st+dt][j]);
				zt[j]=(z[st][j]-z[st+dt][j])*(z[st][j]-z[st+dt][j]);
				s[j][1]=s[j][1]+xt[j];
				s[j][2]=s[j][2]+yt[j];
				s[j][3]=s[j][3]+zt[j];
				s[j][0]=s[j][0]+xt[j]+yt[j]+zt[j];
				// s[j][0]=s[j][0]+xt[j]+yt[j];
				// printf("%lf\n", s[j][0]);
				// printf(" %lf %lf %lf  \n", xt[j], yt[j], zt[j]);
				// printf("Dx %lf %lf %lf\n", x[st][j], x[st+dt][j], x[st][j]-x[st+dt][j]);
				// printf("Dy %lf\n", y[st][j]-y[st+dt][j]);
				// printf("Dz %lf\n", z[st][j]-z[st+dt][j]);
				
				// exit(0);
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
			for (j=0;j<natoms;j++)
			{
				//if (j< N_MOL/2) s1=s1+s[j][dem]/i;
				//if (j>=N_MOL/2) s2=s2+s[j][dem]/i;
				s2=s2+s[j][dem]/i;
			}
			//s1=s1/N_MOL*2.0;
			s2=s2/natoms;
			fprintf(fmsd, "%-15.4lf ", s2);
			//if (dem==0) fprintf(fmsd, "%-15.4lf ", s2/6.0/(0.2*dt)*10); //unit 10-9 m2/s
			if (dem==0) fprintf(fmsd, "%-15.4lf ", s2/6.0/(timestep*dt)*10); //unit 10-9 m2/s
			// if (dem==0) fprintf(fmsd, "%-15.4lf ", s2/4.0/(timestep*dt)*10); //unit 10-9 m2/s
			//printf("%-15.4lf%-15.4lf\n", s1, s2);	
		}
		fprintf(fmsd, "\n");
		//break;
	}
	// printf("Diffusion = %lf\n", s2/4.0/(timestep*dt)*10);
	fclose(fmsd);
	return 0;
}



