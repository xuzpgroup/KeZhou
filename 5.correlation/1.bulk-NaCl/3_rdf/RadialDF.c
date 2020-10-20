#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define MAXNATOM 1000
#define DIM 3
#define Heff 0.00


double 	r1[MAXNATOM][3];
double 	r2[MAXNATOM][3];
double  dis[MAXNATOM][MAXNATOM];
double  box[3];
double  his[1000][2];
double  his0[1000][2];

double DABS(double x)
{
	if (x<0) return -x;
	return x;
}
double Distant_PBC_2D(double *a, double *b, double *cell)
{
    double r;
    double d1, d2;
    d1=DABS(a[0]-b[0]); while (d1>cell[0]) d1=d1-cell[0];
    if (d1>cell[0]/2) d1=cell[0]-d1;
    d2=DABS(a[1]-b[1]); while (d2>cell[1]) d2=d2-cell[1];
    if (d2>cell[1]/2) d2=cell[1]-d2;
    r=sqrt(d1*d1+d2*d2);  
    return r;   
}

double Distant_PBC(double *a, double *b, double *cell)
{
    double r;
    double d1, d2, d3;
    d1=DABS(a[0]-b[0]);	while (d1>cell[0]) d1=d1-cell[0];
    if (d1>cell[0]/2) d1=cell[0]-d1;

    d2=DABS(a[1]-b[1]);	while (d2>cell[1]) d2=d2-cell[1];
    if (d2>cell[1]/2) d2=cell[1]-d2;
    
    d3=DABS(a[2]-b[2]);	while (d3>cell[2]) d3=d3-cell[2];
    if (d3>cell[2]/2) d3=cell[2]-d3;

    r=sqrt(d1*d1+d2*d2+d3*d3);  
    return r;	
}


int main()
{
	int len;
	int step, natoms, n1, n2;
	int i, j;
	double selected;
	int id, itype, imol;
	double chg, vx, vy, vz, xu, yu, zu;
	double xlo, xhi, ylo, yhi, zlo, zhi;

	double x, y, z;
	char header_char[300];
	char filename[300];
	double rm, rst, delt;
	int	bins;


	FILE *fin;
	FILE *fout;

	printf("Start .....\n");

	box[0]=6.3000000000000000e+01;
	box[1]=1.1154400000000000e+02;
	box[2]=10.2;

	rm=10.0;
	rst=0.0;
	bins=200;
	delt=(rm-rst)/bins;

	for (i=0;i<bins;i++)
	{
		his[i][0]=0.0;
		his[i][1]=0.0;	
	}

	selected=0.2;

	len=0;

	fin=fopen("../ions.dat", "r");
	// fin=fopen("../water.v.xyz", "r");
	fout=fopen("rdf.self.dat", "w");	

	double rho, dv, r0;
	int	   rr;

  	printf("Starting ...\n");	
	while (feof(fin)==0)
	{
		len++;
		if (len>40000) break;
		if (len % 100==0) printf("Now is %d\n\n", len);
		// input water and ions from FILE: ions.dat
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &step);
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &natoms);
		fgets(header_char, 300, fin);
		fscanf(fin, "%lf%lf\n", &xlo, &xhi);
		fscanf(fin, "%lf%lf\n", &ylo, &yhi);
		fscanf(fin, "%lf%lf\n", &zlo, &zhi);
		fgets(header_char, 300, fin);
		box[0]=xhi-xlo;
		box[1]=yhi-ylo;
		box[2]=zhi-zlo;
		
		// natoms=natoms/2;

		n1=0;
		n2=0;
		for (i=0;i<natoms;i++)
		{
			// fscanf(fin, "  %d%lf%lf%lf\n", &tp, &x, &y, &z);
			//  ITEM: ATOMS id type xu yu zu
			fscanf(fin, "%d%d%d%lf", &id, &imol, &itype, &chg);
			fscanf(fin, "%lf%lf%lf", &vx, &vy, &vz);
			fscanf(fin, "%lf%lf%lf\n", &xu, &yu, &zu);			// printf("%d\n", i);
			if (itype==9) // observed/center atom 
			{
				r1[n1][0]=xu;
				r1[n1][1]=yu;
				r1[n1][2]=zu;
				n1++;					
			}
			if (itype==10)
			{
				r2[n2][0]=xu;
				r2[n2][1]=yu;
				r2[n2][2]=zu;
				n2++;					
			}
				
		}
		for (i=0;i<n1;i++)
		{
			for (j=0;j<n2;j++)
			{
				dis[i][j]=Distant_PBC(r1[i], r2[j], box);
				// dis[j][i]=dis[i][j];
			}
		}

		if (DIM ==3) rho=n2/box[0]/box[1]/(box[2]-Heff);
		if (DIM ==2) rho=n2/box[0]/box[1];


		for (i=0;i<bins;i++)
		{
			his0[i][0]=0.0;
			his0[i][1]=0.0;	
		}

		for (i=0;i<n1;i++)
		{
			for (j=0;j<n2;j++)
			{
				//if (i==j) continue;
				r0=dis[i][j]-rst;
				rr=r0/delt;
				if (rr>=bins) continue;
				if (DIM ==3) dv=4.0/3.0*3.1415926*((r0+delt)*(r0+delt)*(r0+delt)-r0*r0*r0);
				if (DIM ==2) dv=3.1415926*((r0+delt)*(r0+delt)-r0*r0);
				his0[rr][0]+=1.0/dv/rho;
				his0[rr][1]+=1.0;
				if (dis[i][j]<0.5)
				{
					printf("Find ... %d %d %lf \n",i ,j, dis[i][j]);
					printf("R1 %lf %lf %lf \n", r1[i][0], r1[i][1], r1[i][2]);
					printf("R2 %lf %lf %lf \n", r2[j][0], r2[j][1], r2[j][2]);
					printf("dis %lf \n", dis[i][j]);
				} 
			}
		}

		for (i=0;i<bins;i++)
		{
			his[i][0]+=(his0[i][0]/n1);
			his[i][1]+=(his0[i][1]/n1);
		}		

	}	
	for (i=0;i<bins;i++)
	{
		his[i][0]/=len;
		his[i][1]/=len;
	}

	for (i=1;i<bins;i++)
	{
		his[i][1]=his[i][1]+his[i-1][1];
	}

	for (i=0;i<bins;i++)
	{
		fprintf(fout, "%12.6lf ", rst+(i+0.5)*delt);
		fprintf(fout, "%12.6lf ", his[i][0]);
		fprintf(fout, "%12.6lf ", his[i][1]);
		fprintf(fout, "%12.6lf ", -log(his[i][0]));
		fprintf(fout, "\n");
	}

	fclose(fout);
	fclose(fin);
	
	printf("totial len=%d\n", len);

	FILE *finfo;
	int p1p;

	finfo=fopen("rdf.info.dat", "w");
	//1 st Peak
	for (i=5;i<bins-1;i++)
	{
		if (his[i][0]>his[i-1][0] && his[i][0]>his[i+1][0] && his[i][0]>his[i-2][0] && his[i][0]>his[i+2][0])
		{
			printf("1st Peak = %lf goo = %lf \n", rst+(i+0.5)*delt, his[i][0]);
			fprintf(finfo, "%lf %lf %lf\n", rst+(i+0.5)*delt, his[i][0], his[i][1]);
			p1p=i;
			break;
		}
	}	
	// 1 st Valley
	for (i=5;i<bins-1;i++)
	{
		if (his[i][0]<his[i-1][0] && his[i][0]<his[i+1][0])
		{
			printf("1st Valley = %lf goo = %lf \n", rst+(i+0.5)*delt, his[i][0]);
			fprintf(finfo, "%lf %lf\n", rst+(i+0.5)*delt, his[i][0]);
			break;
		}
	}
	//2 st Peak
	for (i=p1p+2;i<bins-3;i++)
	{
		if (his[i][0]>his[i-1][0] && his[i][0]>his[i+1][0] && his[i][0]>his[i-2][0] && his[i][0]>his[i+2][0])
		{
			printf("2nd Peak = %lf goo = %lf \n", rst+(i+0.5)*delt, his[i][0]);
			fprintf(finfo, "%lf %lf\n", rst+(i+0.5)*delt, his[i][0]);
			break;
		}
	}
	fclose(finfo);	


	
	return 0;
}



