// Culculate the height of metal ions 

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

double Min(double x, double y)
{
	if (x<y) return x;
	return y;
}

double Distant(double *a, double *b)
{
    double r;
    double d1, d2, d3;
    d1=(a[0]-b[0]);
    d2=(a[1]-b[1]);
    d3=(a[2]-b[2]);
    r=sqrt(d1*d1+d2*d2+d3*d3);  
    return r;
}


int main()
{
	char header_char[300];
	double lattice[3][3];
	double fac;
	int natoms[10], natoms_all;
	int i,j;
	double d1, ave_h, min_d, min_z;

	double r[100][3];
	double ff[3];
	char   flag1, flag2, flag3;

	double ref[3];
	double trans[3];

	// relection and translation operation
	ref[0]=1.0;
	ref[1]=1.0;
	ref[2]=1.0;
	trans[0]=1.23436;
	trans[1]=-0.715958;
	trans[2]=0.0;

	FILE *fin;
	FILE *fout;

	ave_h=0.0;
	min_d=10.0;
	min_z=10.0;

	fin=fopen("./CONTCAR", "r");

	// fout=fopen("ions.tar.dat", "w");
	
	printf("Load ...\n");
	fgets(header_char, 300, fin);
	fscanf(fin, "%lf\n", &fac);
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
		{
			fscanf(fin, "%lf", &lattice[i][j]);
		}
		fscanf(fin, "\n");
	}
	fgets(header_char, 300, fin);
	fscanf(fin, "%d %d\n", &natoms[0], &natoms[1]);
	// natoms[0]=32;
	fgets(header_char, 300, fin);
	fgets(header_char, 300, fin);
	// fgets(header_char, 300, fin);
	natoms_all=natoms[0]+natoms[1];

	printf("Contains %d C ...\n", natoms[0]);
	printf("Total atoms is %d ...\n", natoms_all);

	for (i=0;i<natoms_all;i++)
	{
		for (j=0;j<3;j++)
		{
			fscanf(fin, "%lf", &ff[j]);
		}
		fscanf(fin, "   %c   %c   %c\n", &flag1, &flag2, &flag3);
		// printf("%d %lf %lf %lf %c %c %c\n", i+1, ff[0], ff[1], ff[2], flag1, flag2, flag3);
		// for Direct 
		for (j=0;j<3;j++)
		{
			r[i][j]=lattice[0][j]*ff[0]+lattice[1][j]*ff[1]+lattice[2][j]*ff[2];
		}
		// printf("%d %lf %lf %lf \n", i+1, r[i][0], r[i][1], r[i][2]);

		if (i<natoms[0])
		{
			ave_h+=r[i][2];
		}
	}
	fclose(fin);


	ave_h/=natoms[0];
	printf("The average height of graphen wall is %lf \n", ave_h);
	printf("The postion of M is %lf\n", r[natoms_all-1][2]);
	
	for (i=0;i<natoms[0];i++)
	{
		d1=Distant(r[i], r[natoms[0]]);
		min_d=Min(min_d, d1);
		min_z=Min(min_z, r[i][2]);
	}
	printf(" --- The d_MC(ave) is %lf \n", r[natoms[0]][2]-ave_h);
	printf(" --- The min d_MC is %lf \n", min_d);
	printf("The min_z = %lf, The min z_MC is %lf \n", min_z, r[natoms[0]][2]-min_z);


	


	return 0;
}



