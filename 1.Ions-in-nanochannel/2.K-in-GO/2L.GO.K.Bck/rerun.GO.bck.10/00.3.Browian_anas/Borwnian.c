#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define N_MOL 1
#define MAX_T 5000010
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

void Write_Lt(int intervel, int len, double Totaltime)
{
	char filename[100];
	FILE *ft;
	int i;
	double dr[3];

	sprintf(filename, "Lt.%.0lf.dat", Totaltime);
	ft=fopen(filename, "w");
	for (i=0;i+intervel<len;i++)
	{
		dr[0]=(x[i][0]-x[i+intervel][0]);
		dr[1]=(y[i][0]-y[i+intervel][0]);
		dr[2]=(z[i][0]-z[i+intervel][0]);
		fprintf(ft, "%lf %lf %lf\n", dr[0], dr[1], dr[2]);
	}
	fclose(ft);
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
	double timestep, Totaltime;
	int intervel;
	int id, mol, type;
	double  q, x1, y1, z1, xu, yu, zu, vx, vy, vz; 

	FILE *dumpxu;

	printf("Start .....\n");

	len=0;
	dumpxu=fopen("../ions.dat", "r");
	//dumpxu=fopen("mol.position.xyz", "r");

	//printf("%ld\n", len);
	//fgets(header_char, 300, dumpxu);
	//fgets(header_char, 300, dumpxu);
	//fgets(header_char, 300, dumpxu);
	
	tt=0;
	timestep=0.01; // ## ps
	Totaltime=10; //ps
	intervel=(int)Totaltime/timestep;

	while (feof(dumpxu)==0)
	{
		len++;
		tt++;
		// printf("%ld\n", len);
		for (i=0;i<9;i++)
		{
			fgets(header_char, 300, dumpxu);
		}
	
		for (i=0;i<1;i++)
		{
			fscanf(dumpxu, "%d%d%d",&id, &mol, &type);
			fscanf(dumpxu, "%lf%lf%lf", &xu, &yu, &zu);
			fscanf(dumpxu, "%lf%lf%lf", &vx, &vz, &vz);
			fscanf(dumpxu, "\n");				
			if (type==9)
			{
				x[len-1][0]=xu;
				y[len-1][0]=yu;
				z[len-1][0]=zu;
			}
		}

		if (tt<000) len--;
	}	
	// len--;
	fclose(dumpxu);
	
	printf("Calaltion .....\n");
	printf("Total length is %d\n", len);

	Write_Lt(intervel, len, Totaltime);
	Write_Lt(intervel*2, len, Totaltime*2);
	Write_Lt(intervel*3, len, Totaltime*3);
	Write_Lt(intervel*4, len, Totaltime*4);
	Write_Lt(intervel*5, len, Totaltime*5);
	Write_Lt(intervel*6, len, Totaltime*6);
	Write_Lt(intervel*10, len, Totaltime*10);

	




	return 0;
}



