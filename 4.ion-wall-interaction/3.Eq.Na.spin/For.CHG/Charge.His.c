#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define SQRT3	  1.732050808
#define MAXBINS	  1000
#define MAXATOMS  1000

#define Bohr2A    1.889725992
#define e_units   0.036749309 

int    bins;
double hst, hed, delt;
double his[MAXBINS];

double DABS(double x)
{
	if (x<0) return -x;
	return x;
}
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

int equal(double x, double y)
{
	if (DABS(x-y)<0.01) return 1;
	return 0;
}

void CoordOperation(double lattice[3][3], int *rp, double *pp)
{
	int i, j;
	for (i=0;i<3;i++)
	{
		pp[i]=0;
		for (j=0;j<3;j++)
		{
			pp[i]+=lattice[j][i]*rp[j];
		}
	}
}

// void HistogramIt(double samples1, double samples2, double factor)
// {
// 	int kk1, kk2; 
// 	kk1=(samples1-hst[0])/delt[0];
// 	kk2=(samples2-hst[1])/delt[1];
// 	// printf("%d %d \n", kk1, kk2);
// 	if (equal(samples1, hed[0])) kk1--;
// 	if (equal(samples2, hed[1])) kk2--;
	
	
// 	// printf("%lf %d  jj=%d \n", samples, kk, jj);
// 	if (kk1<bins[0] && kk1>=0 && kk2<bins[1] && kk2>=0)
// 	{
// 		his[kk1][kk2]+=1.0*factor;
// 		// other 
// 	}
// 	else {
// 		printf("---%lf %lf %d %d \n", samples1, samples2, kk1, kk2);
// 		printf("Error when histogram! \n");
// 	}
// }

int main()
{
	int 	i, j, k;
	char 	header_char[300];
	char 	filename[300];
	double 	s1, s2, s3;
	double  dlac[3][3];
	double  R[MAXATOMS][3], p[3];
	int     Tp[MAXATOMS];
	int		NGX, NGY, NGZ, natoms;
	int     rp[3];
	double  total_s, grho, scale, lx, lz, transfer_e;

	FILE 	*fin, *fout;





	// get mass center
	printf(" Load cubic file \n");
	fin=fopen("./CHGCAR.cube", "r");
	fgets(header_char, 300, fin);
	fgets(header_char, 300, fin);
	//
	fscanf(fin, "%d%lf%lf%lf\n", &natoms, &s1, &s2, &s3);
	//
	fscanf(fin, "%d%lf%lf%lf\n", &NGX, &dlac[0][0], &dlac[0][1], &dlac[0][2]);
	fscanf(fin, "%d%lf%lf%lf\n", &NGY, &dlac[1][0], &dlac[1][1], &dlac[1][2]);
	fscanf(fin, "%d%lf%lf%lf\n", &NGZ, &dlac[2][0], &dlac[2][1], &dlac[2][2]);
	//
	bins=NGZ;
	for (i=0;i<NGZ;i++) his[i]=0.0;

	for (i=0;i<natoms;i++)
	{
		fscanf(fin, "%d", &Tp[i]);
		fscanf(fin, "%lf", &s1);
		fscanf(fin, "%lf%lf%lf\n", &R[i][0], &R[i][1], &R[i][2]);	
	}

	total_s=0.0;

	for (i=0;i<NGX;i++)
	{
		for (j=0;j<NGY;j++)
		{
			for (k=0;k<NGZ;k++)
			{	
				rp[0]=i; 
				rp[1]=j; 
				rp[2]=k; 				
				CoordOperation(dlac, rp, p);
				fscanf(fin, "%lf", &grho);
				total_s+=grho;
				his[k]+=grho;
			}
		}
	}

	printf(" Total = %.20lf\n", total_s);

	fclose(fin);


//  output 1D his
	total_s=0.0;
	transfer_e=0.0;
	lx=NGX*dlac[0][0]/Bohr2A;
	lz=NGZ*dlac[2][2]/Bohr2A;
	printf(" ... cell info: lx = %lf, lz = %lf\n", lx, lz);
	scale=1.0/NGX/NGY/NGZ/e_units*(lx*lx*SQRT3/2.0*lz);
	fout=fopen("His.in.z.dat", "w");
	fprintf(fout, "####### %d\n  ####  units e/A\n", bins);
	for (i=0;i<bins;i++)
	{
		fprintf(fout, "%lf  %lf\n", dlac[2][2]*i/Bohr2A, his[i]*scale/(lz/NGZ));
		total_s+=his[i]*scale;
		if (dlac[2][2]*i/Bohr2A>11.60 ) transfer_e+=his[i]*scale;
	}
	fclose(fout);
	printf("Check: Total = %.20lf\n\n", total_s);
	printf("    transfer charge = %.20lf\n", transfer_e);

	return 0;
}



