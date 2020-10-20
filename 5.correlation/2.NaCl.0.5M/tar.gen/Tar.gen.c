#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_ATOMS 400
#define MAX_LEN   10000

#define SQRT3	  1.732050808

#define BINS 32
#define NDUMP 250

int flag[MAX_ATOMS];
int n_flag;
// double rs[MAX_LEN][10][3];


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

void Alloc_Flag()
{
	int i;
	for (i=0;i<MAX_ATOMS;i++)
		flag[i]=0;
	n_flag=2;
	flag[15-1]=1;
	flag[26-1]=2;
	
}
int main()
{
	int 	i, j, len;
	char 	header_char[300];
	char 	filename[300], tp;
	int		natoms, step;
	double 	data[20];
	double  Rlist[MAX_ATOMS][3];
	double  Orgin[3], nr[3], rr[3];
	double  s1, s2;
	FILE 	*fin, *fout;
	int		NOUT;
	double  MaxL[2], MinL[2], MaxLength;
	double 	MC[2];
	int     TotalL;

	NOUT=1000;
	TotalL=1000;
	MaxL[0]=-100;
	MaxL[1]=-100;
	MinL[0]=100;
	MinL[1]=100;
	MC[0]=0.0;
	MC[1]=0.0;

	Alloc_Flag();

	// printf(" Calculate mass center \n");
	fin=fopen("../ions.dat", "r");
	fout=fopen("Tar.dat", "w");

	len=0;
	while (!feof(fin))
	{
		len++;
		if (len % 1000 ==1) printf("Now len = %d\n", len);
		// fscanf(fin, "%d\n", &natoms);
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &step);
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &natoms);
		// printf("Natom = %d \n", natoms);
		for (i=0;i<5;i++) fgets(header_char, 300, fin);

		fprintf(fout, " %d %lf ", len*NDUMP, len*NDUMP*0.001) ;
		for (i=0;i<natoms;i++)
		{
			for (j=0;j<10;j++) fscanf(fin, "%lf", &data[j]);
			fscanf(fin, "\n");
			if (flag[i])
			{
				// fprintf(fout, " %d %lf ", len*NDUMP, len*NDUMP*0.001) 
				// printf("%lf %lf %lf \n", data[7], data[8], data[9]);
				fprintf(fout, "%lf %lf %lf   ", data[7], data[8], data[9]);
			}
		}
		fprintf(fout, "\n");

		if (len*NDUMP*0.001>100.0) break;
	}

	fclose(fin);
	fclose(fout);

	// MC[0]=MC[0]/TotalL;
	// MC[1]=MC[1]/TotalL;

	// printf(" Output tarjectory \n");
	// // output tarjectory
	// fin=fopen("../ions.dat", "r");
	// fout=fopen("Tar.dat", "w");
	// len=0;
	// while (!feof(fin))
	// {
	// 	len++;
	// 	if (len % 100000 ==0) printf("Now len = %d\n", len);
	// 	// fscanf(fin, "%d\n", &natoms);
	// 	for (i=0;i<9;i++)
	// 	{
	// 		fgets(header_char, 300, fin);
	// 	}
	// 	for (i=0;i<9;i++) fscanf(fin, "%lf", &data[i]);
	// 	fscanf(fin, "\n");
	// 	RIons[0]=data[3];
	// 	RIons[1]=data[4];
	// 	RIons[2]=data[5];
	// 	if (len%(NOUT/NDUMP)!=0) continue;
	// 	MaxL[0]=MAX_D(RIons[0], MaxL[0]);
	// 	MaxL[1]=MAX_D(RIons[1], MaxL[1]);
	// 	MinL[0]=MIN_D(RIons[0], MinL[0]);
	// 	MinL[1]=MIN_D(RIons[1], MinL[1]);	
		
		
	// 	fprintf(fout, " %d %lf %lf \n", len*NDUMP, RIons[0]-MC[0], RIons[1]-MC[1]);
	// 	if (len/(NOUT/NDUMP)>TotalL) break;
	// }
	// fclose(fin);


	// MaxLength=MAX_D(MaxL[0]-MinL[0], MaxL[1]-MinL[1]);

	// printf("Suggestion:\n");
	// printf("            The Stx = %lf nm \n", MinL[0]/10.0);
	// printf("            The Sty = %lf nm \n", MinL[1]/10.0);
	// printf("            The L   = %lf nm \n", MaxLength/10.0);
	
	
		
	// fclose(fout);

	

	return 0;
}



