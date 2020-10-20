#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define h 10.0
#define zc 4.697650
#define XC -15.760400   
#define YC 42.000800


int main()
{
	int i, natoms, tp;
	char header_char[300];
	double box[3];
	double x, y, z;
	double dx, dy;



	FILE *fin;
	FILE *fout;

	printf("Start .....\n");

	fin=fopen("./graphene.only.xyz", "r");
	fout=fopen("catom.xyz", "w");
	
	fscanf(fin, "%d\n", &natoms);
	fscanf(fin, "%lf%lf%lf\n", &box[0], &box[1], &box[2]);
	fprintf(fout, "%d\n", natoms*2);
	fprintf(fout, "%lf %lf %lf \n", box[0], box[1], box[2]);
	dx=XC-box[0]*0.7;
	dy=YC-box[1]*0.1;

	for (i=0;i<natoms;i++)
	{
		fscanf(fin, "%d", &tp);
		fscanf(fin, "%lf%lf%lf", &x, &y, &z);
		fscanf(fin, "\n");
		if (zc>0)
		{
			fprintf(fout, "C   %lf %lf %lf\n", x+dx, y+dy, zc);
			fprintf(fout, "C   %lf %lf %lf\n", x+dx, y+dy, zc-h);
		}
		else 
		{
			fprintf(fout, "C   %lf %lf %lf\n", x+dx, y+dy, zc+h);
			fprintf(fout, "C   %lf %lf %lf\n", x+dx, y+dy, zc);			
		}

	}
	
	fclose(fin);
	fclose(fout);
	return 0;
}



