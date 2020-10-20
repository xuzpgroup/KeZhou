#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define N_MOL 1
#define MAX_T 50000
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
	int id, mol, type;
	double  q, x1, y1, z1, xu, yu, zu; 

	FILE *dumpxu;
	FILE *fmsd;
	FILE *fdp, *fz;
	FILE *ftar;
	printf("Start .....\n");

	len=0;
	dumpxu=fopen("../ions.dat", "r");
	//dumpxu=fopen("mol.position.xyz", "r");
	fmsd=fopen("msd.data", "w");
	fdp=fopen("gas.dump", "w");
	fz=fopen("gas.z.dat", "w");
	ftar=fopen("gas.tar.dat", "w");
	//printf("%ld\n", len);
	//fgets(header_char, 300, dumpxu);
	//fgets(header_char, 300, dumpxu);
	//fgets(header_char, 300, dumpxu);
	
	tt=0;
	timestep=0.05; // ## ps
	
	while (feof(dumpxu)==0)
	{
		len++;
		tt++;
		// printf("%ld\n", len);
		for (i=0;i<9;i++)
		{
			fgets(header_char, 300, dumpxu);
		}
		
		fprintf(fdp, "%d\n", 3*N_MOL);
		fprintf(fdp, "step: %d\n", len*1000);

		for (i=0;i<1;i++)
		{
			fscanf(dumpxu, "%d%d%d",&id, &mol, &type);
			fscanf(dumpxu, "%lf%lf%lf%lf", &q, &x1, &y1, &z1);
			fscanf(dumpxu, "%lf%lf%lf", &xu, &yu, &zu);
			fscanf(dumpxu, "\n");				
			if (type==9)
			{
				x[len-1][0]=xu;
				y[len-1][0]=yu;
				z[len-1][0]=zu;
			}
			if (type ==9) fprintf(fdp, "Na ");
			else fprintf(fdp, "O ");
			fprintf(fdp, "%lf %lf %lf\n", xu, yu, zu);
			if (type ==9) fprintf(fz, "%12.8lf %lf\n", len*timestep, zu);
			if (len % 200 ==0) 
			fprintf(ftar, "%lf %lf %lf %lf \n", len*timestep, xu, yu, zu);
		}

		if (tt<000) len--;
	}	
	len--;
	fclose(dumpxu);
	printf("Calaltion .....\n");
	printf("Total length is %d\n", len);
	for (dt=1;dt<=len/2;dt=dt+1)
	{
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
		if ( dt%100==1 )printf("Now is %-12d\n", dt);

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
				//s[j][0]=s[j][0]+xt[j]+yt[j]+zt[j];
				s[j][0]=s[j][0]+xt[j]+yt[j]+zt[j];
				// s[j][0]=s[j][0]+xt[j]+yt[j];
				// printf("%lf\n", s[j][0]);
				// printf(" %lf %lf %lf  \n", xt[j], yt[j], zt[j]);
				// printf("Dx %lf %lf %lf\n", x[st][j], x[st+dt][j], x[st][j]-x[st+dt][j]);
				// printf("Dy %lf\n", y[st][j]-y[st+dt][j]);
				// printf("Dz %lf\n", z[st][j]-z[st+dt][j]);
				
				// exit(0);
				// printf("st=%d dt=%d %lf\n", st, dt, xt[j]+yt[j]+zt[j]);
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
			//printf("%-15.4lf%-15.4lf\n", s1, s2);	
		}
		fprintf(fmsd, "\n");
		// break;
	}
	
	fclose(fmsd);
	fclose(fdp);
	return 0;
}



