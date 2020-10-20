#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

// #define RTCOUTPUT
// #define RTCCAL
#define MAXNATOM 10000
#define MAXWT 2000
#define Rcut 3.465
#define Rcut2nd 5.000
#define HisItems 12
#define MAXBINS 500
#define MAXT 500000
#define NPLANE 3312 //for GO
// #define NPLANE 2760 //for G

#define ELEOUT 0


double 	r[MAXNATOM][3];
int 	itype[MAXNATOM];
double  dis[MAXNATOM];
int     disID[MAXNATOM];
int	    shell_list[16];
 		// ## 0 -- num
 		// 1,2,3.. ID



double DABS(double x)
{
	if (x<0) return -x;
	return x;
}

double DMAX(double x, double y)
{
	if (x>y) return x;
	return y;
}

double pbc(double ddx, double lx)
{
	double dx;
	dx=ddx;
	while (dx>0.5*lx) dx=dx-lx;
	while (dx<-0.5*lx) dx=dx+lx;
	return dx;
}
void PutIntoBox(double *r, double *br, double *box)
{
	br[0]=pbc(r[0], box[0]);
	br[1]=pbc(r[1], box[1]);
	br[2]=pbc(r[2], box[2]);
}

void Vec_Minus(double *a, double *b, double *s, int normlized, double *box)
{
	double l;
	// s=a-b
	s[0]=pbc(a[0]-b[0], box[0]);
	s[1]=pbc(a[1]-b[1], box[1]);
	s[2]=pbc(a[2]-b[2], box[2]);	
	l=sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
	if (normlized)
	{
		s[0]/=l;
		s[1]/=l;
		s[2]/=l;
	}
}

void Vec_Plus(double *a, double *b, double *s, int normlized)
{
	double l;
	// s=a+b
	s[0]=a[0]+b[0];
	s[1]=a[1]+b[1];
	s[2]=a[2]+b[2];	
	l=sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
	if (normlized)
	{
		s[0]/=l;
		s[1]/=l;
		s[2]/=l;
	}
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
	int step, natoms, nwater;
	int i, j, st, tp;
	double selected;
	int id, imol;
	double box[3];
	// for life time 
	double timestep;
	int shell_life;
	int new_shell;
	double slife_time, slife_OH;
	int nconst, nconst_OH; // total contruction times 
	int OHshell;
	// for theta 

	int lw[100];
	int lf[100];

	double rrin[3];
	int idwater, jdwater;
	double ss[HisItems];

	char header_char[300];
	char filename[300];
	double data[8];
	int gid, gmol, gtype;
	double gq, gr[3], gru[3];

	double sRp;			// rattling of ion within its solvation shell
						// the distortion response of the solvation shell to the external perturbation 
						// P=<|R_ion-R_com|>
	double Rcom[3];	

	FILE *fin;
	FILE *flw;
	FILE *flf;


	box[0]=63;
	box[1]=111.544;
	box[2]=10.0;

	// important paramters
	timestep=0.050; // ps
	
	printf("Start .....\n");
	


	len=0;
	//fin=fopen("../see.xyz", "r");
	fin=fopen("../sys.unwrap.xyz", "r");
	// fions=fopen("../ions.dat", "r");

	flw=fopen("./list.water.dat", "w");
	flf=fopen("./list.function.dat", "w");
	
	

	while (feof(fin)==0)
	{
		// if (len>40000) break;
		len++;
		if (len % 1000 ==0) printf("Now is %d\n", step);
	
		// input water and ions from FILE: sys.unwrap.xyz
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &step);
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &natoms);
		for (i=0;i<5;i++) fgets(header_char, 300, fin);
		
		nwater=0; // including to O atom in -OH
		for (i=0;i<natoms;i++)
		{
			// ITEM: ATOMS id type xu yu zu
			fscanf(fin, "%d%d", &id, &itype[i]);
			fscanf(fin, "%lf%lf%lf\n", &r[i][0], &r[i][1], &r[i][2]);
			if (itype[i]==8 || itype[i]==5 || itype[i]==3)
			{
				nwater++;
				disID[nwater-1]=i;
			}
			if (itype[i]==9)
			{
				gru[0]=r[i][0];
				gru[1]=r[i][1];
				gru[2]=r[i][2];					
			}	
		}
		// calculate dist to mental ions 
		// fprintf(frcfinfo, "#### %d\n", len);
		
		lw[0]=0;
		lf[0]=0;

		for (i=0;i<nwater;i++)
		{
			id=disID[i];
			dis[i]=Distant_PBC(r[id], gru, box);
			// if (len>MAXT) continue;
			if (dis[i]<Rcut)
			{
				if (itype[id]==8) // is O atom in water 
				{
					lw[0]++;
					lw[lw[0]]=id;
				}
				if (itype[id]==3 || itype[id]==5) // is O atom in function group 
				{
					lf[0]++;
					lf[lf[0]]=id;
				}
			}
		}
		// output list 
		// the list for water
		fprintf(flw, "%d %d\n", len, lw[0]);
		for (i=0;i<lw[0];i++)
		{
			fprintf(flw, "%d\n", lw[i+1]);
		}
		// the list for function group
		fprintf(flf, "%d %d\n", len, lf[0]);
		for (i=0;i<lf[0];i++)
		{
			fprintf(flf, "%d\n", lf[i+1]);
		}
	}	
	// fclose(fions);
	fclose(fin);
	fclose(flw);
	fclose(flf);


	return 0;
}



