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
int	    last_shell_list[16];

short	v[MAXT][MAXWT];

// for histogram
int    bins;
double hst[HisItems], hed[HisItems], delt[HisItems];
double his[MAXBINS][HisItems];
	// 	info
	// 	0	theta 1
	// 	1	cos (theta1)
	// 	2 	theta 2
	// 	3	cos (theta2)
	// 	4	Nc
	// 	5	lifetime 
	//	6	cos (thetaz)
	//  7   span x
	//  8	span y
	//  9   span z
	//  10	Rmo
	//  11  RHO

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

void Vec_Cross(double *a, double *b, double *s, int normlized)
{
	double u1, u2, u3, v1, v2, v3;
	double l;
	u1=a[0];	u2=a[1];	u3=a[2];
	v1=b[0];	v2=b[1];	v3=b[2];
	// s=a*b
	s[0]=u2*v3-u3*v2;
	s[1]=u3*v1-u1*v3;
	s[2]=u1*v2-u2*v1;	
	l=sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
	if (normlized)
	{
		s[0]/=l;
		s[1]/=l;
		s[2]/=l;
	}
}

double Vec_Dot(double *a, double *b)
{
	double s;
	s=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	return s;
}


void swapint(int *a, int *b)
{
    int ss;
    ss=*a;     *a=*b;    *b=ss;
}
void swapdouble(double *a, double *b)
{
    double ss;
    ss=*a;     *a=*b;    *b=ss;
}

void swapdata(int i, int j)    
{  
	swapint(&disID[i], &disID[j]);
	swapdouble(&dis[i], &dis[j]);
    return ;  
}    

void quicksort(int begin, int end)  
{  
    int i, j;  
  
    if(begin < end)  
    {  
        i = begin + 1;   
        j = end;        
            
        while(i < j)  
        {  
            // if(rmol[i] > rmol[begin]) 
            if(dis[i] > dis[begin]) 
            {  
                // swap(&array[i], &array[j]);
                swapdata(i, j);
                j--;  
            }  
            else  
            {  
                i++;
            }  
        }  
        // if(rmol[i] >= rmol[begin])   
        if(dis[i] >= dis[begin]) 
        {  
            i--;  
        }  
  
        // swap(&array[begin], &array[i]);
        swapdata(begin, i);  
        quicksort(begin, i);  
        quicksort(j, end);  
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

int Constrct_New_shell()
{
	int i,j;
	int n, temp;
	if (shell_list[0]!=last_shell_list[0]) return 1;
	//sort according to ID
	n=shell_list[0];
	for (i=1;i<=n;i++)
	{
		for (j=i+1;j<=n;j++)
		{
			if (shell_list[i]>shell_list[j])
			{
				temp=shell_list[i];
				shell_list[i]=shell_list[j];
				shell_list[j]=temp;
			}
		}
	}
	
	n=last_shell_list[0];
	for (i=1;i<=n;i++)
	{
		for (j=i+1;j<=n;j++)
		{
			if (last_shell_list[i]>last_shell_list[j])
			{
				temp=last_shell_list[i];
				last_shell_list[i]=last_shell_list[j];
				last_shell_list[j]=temp;
			}
		}
	}

	for (i=0;i<=n;i++)
	{
		if (shell_list[i]!=last_shell_list[i])
			return 1;
	}
	return 0;
}

void HistogramIt(double samples, int jj)
{
	int kk; 

	kk=(samples-hst[jj])/delt[jj];
	if (kk<bins)
	{
		his[kk][jj]+=1.0;
	}
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
	double dir1[3], dir2[3], wp[3], mo[3]; 
	double dx, dy, dz;
	double theta2, theta1, cos1, cos2;
	double thetaz, cosz;
	double drmax[3];
	double avespan[3];

	int n2nd;

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
	FILE *fions;
	FILE *fNc;
	FILE *fshell;
	FILE *fshell2nd;
	FILE *flife;
	FILE *fcos;
	FILE *ftheta;
	FILE *fhisnc;
	FILE *fhislife;
	FILE *fhisthetaz;
	FILE *fspan;
	FILE *frcfinfo;
	FILE *fOH;
	FILE *fRmo;

	box[0]=63;
	box[1]=111.544;
	box[2]=10.0;

	// important paramters
	timestep=0.050; // ps
	
	printf("Start .....\n");

	hst[0]=0.0;	hed[0]=180.0;
	hst[1]=-1;	hed[1]=1.0;
	hst[2]=0.0;	hed[2]=180.0;
	hst[3]=-1;	hed[3]=1.0;
	hst[4]=0;	hed[4]=100;
	hst[5]=0;	hed[5]=100;
	hst[6]=0.0;	hed[6]=90.0;
	hst[7]=0.0;	hed[7]=10.0;
	hst[8]=0.0;	hed[8]=10.0;
	hst[9]=0.0;	hed[9]=10.0;
	hst[10]=1.5;	hed[10]=3.50;
	hst[11]=1.5;	hed[11]=3.50;
	
	
	
	
	bins=100;
	for (i=0;i<HisItems;i++) delt[i]=(hed[i]-hst[i])/bins;
	for (i=0;i<bins;i++)
	{
		for (j=0;j<HisItems;j++)
			his[i][j]=0.0;
	}

	len=0;
	//fin=fopen("../see.xyz", "r");
	fin=fopen("../sys.unwrap.xyz", "r");
	// fions=fopen("../ions.dat", "r");

	fNc=fopen("./info.Nc.dat", "w");
	fshell=fopen("./shell.str.xyz", "w");
	fshell2nd=fopen("./shell.str.2nd.xyz", "w");
	flife=fopen("./shell.life.info.dat", "w");
	fcos=fopen("./his.cos.dat", "w");
	ftheta=fopen("./his.theta.dat", "w");
	fhisnc=fopen("./his.nc.dat", "w");
	fhislife=fopen("./his.life.dat", "w");
	fhisthetaz=fopen("./his.thetaz.dat", "w");
	fspan=fopen("./his.span.dat", "w");
	fRmo=fopen("./his.Rmo.dat", "w");
	frcfinfo=fopen("./RCF.info.dat", "w");
	fOH=fopen("./OH.dat", "w");
	
	
		
	
	last_shell_list[0]=0;

	shell_life=0;
	slife_time=0;
	nconst=0;
	slife_OH=0;
	nconst_OH=0;
	OHshell=0;
	avespan[0]=0.0;	avespan[1]=0.0;	avespan[2]=0.0;
	sRp=0.0;
	
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
			if (itype[i]==8 || itype[i]==5)
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
		for (i=0;i<nwater;i++)
		{
			id=disID[i];
			dis[i]=Distant_PBC(r[id], gru, box);
#ifdef RTCOUTPUT
			for RCT claclulation 
			if (len>MAXT) continue;
			v[len-1][i]=0;
			if (dis[i]<Rcut)
			{
				v[len-1][i]=1;
			}
			fprintf(frcfinfo, "%d %d\n", i+1, v[len-1][i]);
#endif
		}

		// nwater=natoms/3;
		quicksort(0, nwater-1);
		//for (i=0;i<4;i++) printf("%lf   ", dis[i]);
		//printf("\n");
		
		// cal the shell list
		shell_list[0]=0;
		OHshell=0;
		for (i=0;i<20;i++)
		{
			if (dis[i]<Rcut)
			{
				shell_list[0]++;
				shell_list[shell_list[0]]=disID[i];
				if (itype[disID[i]]==5) OHshell=1;
				// histogram Rmo
				if (itype[disID[i]]==8) HistogramIt(dis[i], 10);
				if (itype[disID[i]]==5) HistogramIt(dis[i], 11);	
			}
			if (dis[i]>Rcut)
			{
				break;
			}
		}
		fprintf(fNc, "%12.6lf  %d\n", timestep*(len-1), shell_list[0]);
		HistogramIt(shell_list[0], 4);

		//printf("%d  %d\n", len, shell_list[0]);
		
		if (i>=13) printf("ERROR when construction the shell list\n");
		
		// cal life time, constrct new shell or not?
		new_shell=Constrct_New_shell();
		slife_time++;
		if (OHshell) slife_OH++;
		if (!new_shell)
		{
			shell_life++;
		}
		else 
		{
			// slift_time=+shell_life;
			nconst++;
			if (OHshell) 
			{
				nconst_OH++;
				fprintf(fOH, "%d %d %d\n", len, nconst-1, shell_life);
			}
			// update shell list 
			last_shell_list[0]=shell_list[0];
			for (i=1;i<=shell_list[0];i++)	last_shell_list[i]=shell_list[i];
			// histogram the shell life time 
			HistogramIt(shell_life, 5);
			shell_life=1;
		}
		// output shell stucture 
		fprintf(fshell, "%d\n", shell_list[0]*3+1);
		fprintf(fshell, "time=%d\n", len);
		PutIntoBox(gru, rrin, box);
		// fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 9, rrin[0], rrin[1], rrin[2]);
		fprintf(fshell, " Na %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]);
		
		for (i=1;i<=shell_list[0];i++) 
		{	
			idwater=shell_list[i];
			if (itype[idwater]==8) // is water molecue
			{
				PutIntoBox(r[idwater], rrin, box);
				if (!ELEOUT) fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 8, rrin[0], rrin[1], rrin[2]);
				if ( ELEOUT) fprintf(fshell, " OW %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]);
				PutIntoBox(r[idwater-1], rrin, box);
				if (!ELEOUT) fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 7, rrin[0], rrin[1], rrin[2]);
				if ( ELEOUT) fprintf(fshell, " HW %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]);
				PutIntoBox(r[idwater-2], rrin, box);
				if (!ELEOUT) fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 7, rrin[0], rrin[1], rrin[2]);
				if ( ELEOUT) fprintf(fshell, " HW %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]);
			}
			else if (itype[idwater]==5) // is OH
			{
				PutIntoBox(r[idwater], rrin, box);
				if (!ELEOUT) fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 5, rrin[0], rrin[1], rrin[2]); // OH
				if ( ELEOUT) fprintf(fshell, " OH %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]); // OH
				PutIntoBox(r[idwater-1], rrin, box);
				if (!ELEOUT) fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 4, rrin[0], rrin[1], rrin[2]); // CF
				if ( ELEOUT) fprintf(fshell, " CF %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]); // CF
				PutIntoBox(r[idwater+1], rrin, box);
				if (!ELEOUT) fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 6, rrin[0], rrin[1], rrin[2]); // HO			
				if ( ELEOUT) fprintf(fshell, " HO %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]); // HO			
			}
		}

		// output 2nd shell stucture 
		n2nd=0;
		for (i=0;i<100;i++)
		{
			if (dis[i]<Rcut2nd )
			{
				n2nd++;
			}
			if (dis[i]>Rcut2nd)
			{
				break;
			}
		}		
		fprintf(fshell2nd, "%d\n",n2nd*3+1);
		fprintf(fshell2nd, "time=%d\n", len);
		PutIntoBox(gru, rrin, box);
		// fprintf(fshell2nd, "%4d%12.6lf%12.6lf%12.6lf\n", 9, rrin[0], rrin[1], rrin[2]);
		fprintf(fshell2nd, " Na %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]);
		for (i=0;i<100;i++)
		{
			if (dis[i]<Rcut2nd )
			{
				idwater=disID[i];
				//printf("%d %d\n", idwater, itype[idwater]);
				if (itype[idwater]==8) // is water molecue
				{
					PutIntoBox(r[idwater], rrin, box);
					if (!ELEOUT) fprintf(fshell2nd, "%4d%12.6lf%12.6lf%12.6lf\n", 8, rrin[0], rrin[1], rrin[2]);
					if ( ELEOUT) fprintf(fshell2nd, " OW %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]);
					PutIntoBox(r[idwater-1], rrin, box);
					if (!ELEOUT) fprintf(fshell2nd, "%4d%12.6lf%12.6lf%12.6lf\n", 7, rrin[0], rrin[1], rrin[2]);
					if ( ELEOUT) fprintf(fshell2nd, " HW %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]);
					PutIntoBox(r[idwater-2], rrin, box);
					if (!ELEOUT) fprintf(fshell2nd, "%4d%12.6lf%12.6lf%12.6lf\n", 7, rrin[0], rrin[1], rrin[2]);
					if ( ELEOUT) fprintf(fshell2nd, " HW %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]);
				}
				else if (itype[idwater]==5) // is OH
				{
					PutIntoBox(r[idwater], rrin, box);
					if (!ELEOUT) fprintf(fshell2nd, "%4d%12.6lf%12.6lf%12.6lf\n", 5, rrin[0], rrin[1], rrin[2]); // OH
					if ( ELEOUT) fprintf(fshell2nd, " OF %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]); // OH
					PutIntoBox(r[idwater-1], rrin, box);
					if (!ELEOUT) fprintf(fshell2nd, "%4d%12.6lf%12.6lf%12.6lf\n", 4, rrin[0], rrin[1], rrin[2]); // CF
					if ( ELEOUT) fprintf(fshell2nd, " CF %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]); // CF
					PutIntoBox(r[idwater+1], rrin, box);
					if (!ELEOUT) fprintf(fshell2nd, "%4d%12.6lf%12.6lf%12.6lf\n", 6, rrin[0], rrin[1], rrin[2]); // HO			
					if ( ELEOUT) fprintf(fshell2nd, " HO %12.6lf%12.6lf%12.6lf\n", rrin[0], rrin[1], rrin[2]); // HO			
				}			
			}
			if (dis[i]>Rcut2nd)
			{
				break;
			}
		}		



		// cal the theta 2
		// printf("#################\n");
		for (i=1;i<=shell_list[0];i++)
		{
			idwater=shell_list[i];
			if (itype[idwater]==5) continue; // OH 
			// printf("water O %lf  %lf  %lf  \n", r[idwater][0], r[idwater][1], r[idwater][2]);
			// printf("water H %lf  %lf  %lf  \n", r[idwater-1][0], r[idwater-1][1], r[idwater-1][2]);
			// printf("water H %lf  %lf  %lf  \n", r[idwater-2][0], r[idwater-2][1], r[idwater-2][2]);
			Vec_Minus(r[idwater-1], r[idwater], dir1, 0, box);
			Vec_Minus(r[idwater-2], r[idwater], dir2, 0, box);
			// printf("Dir1  %lf  %lf  %lf  \n", dir1[0], dir1[1], dir1[2]);
			// printf("Dir2  %lf  %lf  %lf  \n", dir2[0], dir2[1], dir2[2]);
			Vec_Plus(dir1, dir2, wp, 1);
			// printf("Iwater=%d\n", idwater);
			// printf("wp %lf  %lf  %lf  \n", wp[0], wp[1], wp[2]);
			Vec_Minus(r[idwater], gru, mo, 1, box);
			// printf("mo %lf  %lf  %lf  \n", mo[0],mo[1], mo[2]);
			cos2=Vec_Dot(wp, mo);
			theta2=acos(cos2)*180.0/3.1415926;
			// printf("theta2 %lf  %lf\n", cos2, theta2);
			cosz=DABS(mo[2]);
			thetaz=acos(cosz)*180.0/3.1415926;
			// printf("thetaz %lf  %lf\n", cosz, thetaz);
			// histogram theta2 and cos		
			HistogramIt(cos2, 3);
			HistogramIt(theta2, 2);	
			// histogram thetaz
			HistogramIt(thetaz, 6);	
		}	
		// cal <|R_ion-R_com|>
		Rcom[0]=0.0;	
		Rcom[1]=0.0;	
		Rcom[2]=0.0;
		for (i=1;i<=shell_list[0];i++)
		{
			idwater=shell_list[i];
			Vec_Minus(r[idwater], gru, mo, 0, box);
			Vec_Plus(mo, Rcom, Rcom, 0);
		}	
		Rcom[0]/=shell_list[0];
		Rcom[1]/=shell_list[0];
		Rcom[2]/=shell_list[0];
		sRp+=sqrt(Rcom[0]*Rcom[0]+Rcom[1]*Rcom[1]+Rcom[2]*Rcom[2]);

		
		// cal the theta 1
		for (i=1;i<=shell_list[0];i++)
		{
			idwater=shell_list[i];
			for (j=i+1;j<=shell_list[0];j++)
			{
				jdwater=shell_list[j];
				Vec_Minus(r[idwater], gru, dir1, 1, box);
				Vec_Minus(r[jdwater], gru, dir2, 1, box);
				cos1=Vec_Dot(dir1, dir2);
				theta1=acos(cos1)*180.0/3.1415926;
				// printf("theta1 %lf  %lf\n", cos1, theta1);
				// histogram theta2 and cos		
				HistogramIt(cos1, 1);
				HistogramIt(theta1, 0);

			}
		}
		// cal space span 
		drmax[0]=0.0;	
		drmax[1]=0.0;	
		drmax[2]=0.0;	
		for (i=1;i<=shell_list[0];i++)
		{
			idwater=shell_list[i];
			Vec_Minus(r[idwater], gru, dir1, 0, box);
			dx=DABS(dir1[0]);
			dy=DABS(dir1[1]);
			dz=DABS(dir1[2]);
			drmax[0]=DMAX(dx, drmax[0]);
			drmax[1]=DMAX(dy, drmax[1]);
			drmax[2]=DMAX(dz, drmax[2]);	
			for (j=i+1;j<=shell_list[0];j++)
			{
				jdwater=shell_list[j];
				Vec_Minus(r[idwater], r[jdwater], dir1, 0, box);
				dx=DABS(dir1[0]);
				dy=DABS(dir1[1]);
				dz=DABS(dir1[2]);
				// printf("dr %lf  %lf  %lf\n", dir1[0], dir1[1], dir1[2]);
				// printf("dr %lf  %lf  %lf\n", dx, dy, dz);
				drmax[0]=DMAX(dx, drmax[0]);
				drmax[1]=DMAX(dy, drmax[1]);
				drmax[2]=DMAX(dz, drmax[2]);	
			}
		}
		// printf("span %lf  %lf  %lf\n", drmax[0], drmax[1], drmax[2]);

		// histogram span		
		HistogramIt(drmax[0], 7);
		HistogramIt(drmax[1], 8);
		HistogramIt(drmax[2], 9);
		avespan[0]+=drmax[0];
		avespan[1]+=drmax[1];
		avespan[2]+=drmax[2];
		
		
		//if (len==700) break;
	}	
	// fclose(fions);
	fclose(fin);
	fclose(fshell);
	fclose(fshell2nd);
	fclose(frcfinfo);
	fclose(fOH);
	
	// normalization 
	// the area is zero 
	for (j=0;j<HisItems;j++) ss[j]=0.0;
	for (i=0;i<bins;i++) 
	{
		for (j=0;j<HisItems;j++)
		{
			// his[i][j]/=len;
			ss[j]+=his[i][j]*delt[j];
		}
	}
	for (i=0;i<bins;i++) 
	{
		for (j=0;j<HisItems;j++)
		{
			his[i][j]/=ss[j];
		}
	}
	// his for cos
	fprintf(fcos, "## bin = %d\n", bins);
	fprintf(fcos, "## cos()  ## pcos(o1) ## pcos(o2) \n");
	for (i=0;i<bins;i++)
	{
		fprintf(fcos, "%12.6lf ", hst[1]+(i+0.5)*delt[1]);
		fprintf(fcos, "%12.6lf ", his[i][1]);
		fprintf(fcos, "%12.6lf ", his[i][3]);
		fprintf(fcos, "\n");
	}
	fclose(fcos);	
	// his for theta
	fprintf(ftheta, "## bin = %d\n", bins);
	fprintf(ftheta, "## theta  ## p(o1) ## p(o2) \n");
	for (i=0;i<bins;i++)
	{
		fprintf(ftheta, "%12.6lf ", hst[0]+(i+0.5)*delt[0]);
		fprintf(ftheta, "%12.6lf ", his[i][0]);
		fprintf(ftheta, "%12.6lf ", his[i][2]);
		fprintf(ftheta, "\n");
	}
	fclose(ftheta);	
	// his for Nc
	fprintf(fhisnc, "## bin = %d\n", bins);
	fprintf(fhisnc, "## Nc  ## p(Nc) \n");
	//for (i=0;i<bins;i++)
	double avenc;
	avenc=0;
	for (i=2;i<12;i++)
	{
		fprintf(fhisnc, "%12.0lf ", hst[4]+(i+0.0)*delt[4]);
		fprintf(fhisnc, "%12.6lf ", his[i][4]);
		fprintf(fhisnc, "\n");
		avenc+=(hst[4]+(i+0.0)*delt[4])*his[i][4];
	}
	fclose(fhisnc);	
	// his for lifetime
	fprintf(fhislife, "## bin = %d\n", bins);
	fprintf(fhislife, "## bins ##l ifetime  ## p(lifetime) \n");
	for (i=0;i<bins;i++)
	{
		fprintf(fhislife, "%4d %12.6lf ", i+1, (hst[5]+(i+1.0)*delt[5])*timestep);
		fprintf(fhislife, "%12.6lf ", his[i][5]);
		fprintf(fhislife, "\n");
	}
	fclose(fhislife);		
	// his for theta
	fprintf(fhisthetaz, "## bin = %d\n", bins);
	fprintf(fhisthetaz, "## thetaz  ## p(oz)\n");
	for (i=0;i<bins;i++)
	{
		fprintf(fhisthetaz, "%12.6lf ", hst[6]+(i+0.5)*delt[6]);
		fprintf(fhisthetaz, "%12.6lf ", his[i][6]);
		fprintf(fhisthetaz, "\n");
	}
	fclose(fhisthetaz);	
	// his for span 
	fprintf(fspan, "## bin = %d\n", bins);
	fprintf(fspan, "## dr  ## p(X) p(Y) p(Z)\n");
	for (i=0;i<bins;i++)
	{
		fprintf(fspan, "%12.6lf ", hst[7]+(i+0.5)*delt[7]);
		fprintf(fspan, "%12.6lf ", his[i][7]);
		fprintf(fspan, "%12.6lf ", his[i][8]);
		fprintf(fspan, "%12.6lf ", his[i][9]);
		
		fprintf(fspan, "\n");
	}
	fclose(fspan);	
	// for Rmo Rho
	fprintf(fRmo, "## bin = %d\n", bins);
	fprintf(fRmo, "## r  ## pMO(r) pHO(r)\n");
	for (i=0;i<bins;i++)
	{
		fprintf(fRmo, "%12.6lf ", hst[10]+(i+0.5)*delt[10]);
		fprintf(fRmo, "%12.6lf ", his[i][10]);
		fprintf(fRmo, "%12.6lf ", his[i][11]);
		
		fprintf(fRmo, "\n");
	}
	fclose(fRmo);		
	


	printf("totial len=%d\n", len);
	printf("totial Step=%d\n", step);
	printf("Averagge lifetime %lf\n", slife_time/nconst*timestep);
	printf("Averagge lifetime for -OH%lf\n", slife_OH/nconst_OH*timestep);
	printf("Output ave info\n");
	fprintf(flife, "%lf\n", slife_time/nconst*timestep);
	fprintf(flife, "%lf\n", slife_OH/nconst_OH*timestep);
	fprintf(flife, "%lf\n", avenc);
	fprintf(flife, "%lf %lf %lf \n", avespan[0]/len, avespan[1]/len, avespan[2]/len);
	fprintf(flife, "%lf\n", sRp/len);
	fprintf(flife, "%lf\n", slife_time);
	fprintf(flife, "%lf\n", slife_OH);
	
	fclose(flife);

#ifdef RTCCAL
	// for calculation Residence correlation function 
	double Maxtime=50.0;
	int max_ser, k, jj;
	double tao, summation, s1;
	FILE *frtd;

	frtd=fopen("./RCF.dat", "w");
	
	max_ser=floor(Maxtime/timestep);
	for (i=0;i<=max_ser;i++)
	{
		tao=i*timestep;
		summation=0.0;
		jj=0;
		for (j=0;j<len;j++)
		{
			if (j+i<len)
			{
				jj++;
				s1=0;
				for (k=0;k<nwater;k++)
				{
					s1+=v[j][k]*v[j+i][k];
					// if (v[j][k]) printf(" ####\n");
					summation+=v[j][k]*v[j+i][k];	
				}
				// printf("%lf\n", s1);
			}
		}
		summation=summation/jj;
		fprintf(frtd, "%lf %lf\n", i*timestep, summation);
		// printf("%lf  %d\n", summation, jj);
	}
	fclose(frtd);
#endif 


	return 0;
}



