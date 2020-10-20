#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define MAXNATOM 20000
#define Rcut 3.3 
#define HisItems 7
#define MAXBINS 500


double 	r[MAXNATOM][3];
double  dis[MAXNATOM];
int     disID[MAXNATOM];
int	    shell_list[16];
 		// ## 0 -- num
 		// 1,2,3.. ID
int	    last_shell_list[16];

// for histogram
int    bins;
double hst[HisItems], hed[HisItems], delt[HisItems];
double his[MAXBINS][HisItems];
	// 	info
	// 	0	theta 1
	// 	1	cos (theta1)
	// 	2 	theta 2
	// 	3	cos (theta2)
	// 	4	lifetime 
	// 	5	Nc
	//	6	cos (thetaz)

double DABS(double x)
{
	if (x<0) return -x;
	return x;
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

void swap(int i, int j)
{
    int temp;
    double tempdou;
    temp=disID[i];		disID[i]=disID[j];		disID[j]=temp;
	tempdou=dis[i];		dis[i]=dis[j];		dis[j]=tempdou;
}

void quicksort(int left, int right)
{
    int pivot, i, j;

    if (left >= right) { return; }
    pivot = dis[left];
    i = left + 1;
    j = right;
    while (1)
    {
        while (i <= right)
        {
            if (dis[i] > pivot)
            {
                break;
            }
            i++;
        }
        while (j > left)
        {
            if (dis[j] < pivot)
            {
                break;
            }
            j--;
        }
        if (i > j) { break; }
        swap(i, j);

    }

    swap(left, j);

    quicksort(left, j - 1);
    quicksort(j + 1, right);
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
	int i, j, st, tp, temp;
	double selected;
	int id, itype, imol;
	double box[3];
	// for life time 
	double timestep;
	int shell_life;
	int new_shell;
	double slift_time;
	int nconst; // total contruction times 
	// for theta 
	double dir1[3], dir2[3], wp[3], mo[3]; 
	double theta2, theta1, cos1, cos2;
	double thetaz, cosz;

	double rrin[3];
	int idwater, jdwater;
	double ss[HisItems];

	char header_char[300];
	char filename[300];
	double data[8];
	int gid, gmol, gtype;
	double gq, gr[3], gru[3];
		

	FILE *fin;
	FILE *fions;
	FILE *fNc;
	FILE *fshell;
	FILE *flife;
	FILE *fcos;
	FILE *ftheta;
	FILE *fhisnc;
	FILE *fhislife;
	FILE *fhisthetaz;

	box[0]=63;
	box[1]=111.544;
	box[2]=6.8;
	timestep=0.250; // ps
	
	printf("Start .....\n");

	hst[0]=0.0;	hed[0]=180.0;
	hst[1]=-1;	hed[1]=1.0;
	hst[2]=0.0;	hed[2]=180.0;
	hst[3]=-1;	hed[3]=1.0;
	hst[4]=0;	hed[4]=100;
	hst[5]=0;	hed[5]=100;
	hst[6]=0.0;	hed[6]=90.0;
	
	
	bins=100;
	for (i=0;i<HisItems;i++) delt[i]=(hed[i]-hst[i])/bins;
	for (i=0;i<bins;i++)
	{
		for (j=0;j<HisItems;j++)
			his[i][j]=0.0;
	}

	len=0;

	fin=fopen("../sys.unwrap.xyz", "r");
	fions=fopen("../ions.dat", "r");
	fNc=fopen("./info.Nc.dat", "w");
	fshell=fopen("./shell.str.xyz", "w");
	flife=fopen("./shell.life.info.dat", "w");

	fcos=fopen("./his.cos.dat", "w");
	ftheta=fopen("./his.theta.dat", "w");
	fhisnc=fopen("./his.nc.dat", "w");
	fhislife=fopen("./his.life.dat", "w");
	fhisthetaz=fopen("./his.thetaz.dat", "w");
	
	
	
	
	last_shell_list[0]=0;

	shell_life=0;
	slift_time=0;
	nconst=0;
	
	while (feof(fin)==0)
	{
		len++;
		if (len % 1000 ==0) printf("Now is %d\n", step);
		// input ions 
		// ITEM: ATOMS id mol type q x y z xu yu zu
		for (i=0;i<9;i++) fgets(header_char, 300, fions);
		fscanf(fions, "%d%d%d",  &id, &imol, &itype);
		fscanf(fions, "%lf", &gq);
		fscanf(fions, "%lf%lf%lf", &gr[0], &gr[1], &gr[0]);
		fscanf(fions, "%lf%lf%lf\n", &gru[0], &gru[1], &gru[2]);
		for (temp=0;temp<10*4;temp++) fgets(header_char, 300, fions);
		// input water 
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &step);
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &natoms);
		for (i=0;i<5;i++) fgets(header_char, 300, fin);
			
		for (i=0;i<natoms;i++)
		{
			// ITEM: ATOMS id type xu yu zu
			fscanf(fin, "%d%d", &id, &itype);
			fscanf(fin, "%lf%lf%lf\n", &r[i][0], &r[i][1], &r[i][2]);
			if (i%3==2)
			{
				disID[i/3]=i;
				dis[i/3]=Distant_PBC(r[i], gru, box);
			}
			
		}
		nwater=natoms/3;
		quicksort(0, nwater-1);
		
		// cal the shell list
		shell_list[0]=0;
		for (i=0;i<20;i++)
		{
			if (dis[i]<Rcut)
			{
				shell_list[0]++;
				shell_list[shell_list[0]]=disID[i];
			}
			if (dis[i]>Rcut)
			{
				break;
			}
		}
		fprintf(fNc, "%12.6lf  %d\n", timestep*(len-1), shell_list[0]);
		HistogramIt(shell_list[0], 4);

		if (i>11) printf("ERROR when construction the shell list, Nc=%d\n", shell_list[0]);
		
		// cal life time, constrct new shell or not?
		new_shell=Constrct_New_shell();
		slift_time++;
		if (!new_shell)
		{
			shell_life++;
		}
		else 
		{
			// slift_time=+shell_life;
			nconst++;
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
		fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 9, rrin[0], rrin[1], rrin[2]);
		for (i=1;i<=shell_list[0];i++) 
		{	
			idwater=shell_list[i];
			PutIntoBox(r[idwater], rrin, box);
			fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 8, rrin[0], rrin[1], rrin[2]);
			PutIntoBox(r[idwater-1], rrin, box);
			fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 7, rrin[0], rrin[1], rrin[2]);
			PutIntoBox(r[idwater-2], rrin, box);
			fprintf(fshell, "%4d%12.6lf%12.6lf%12.6lf\n", 7, rrin[0], rrin[1], rrin[2]);
		}
		// cal the theta 2
		// printf("#################\n");
		for (i=1;i<=shell_list[0];i++)
		{
			idwater=shell_list[i];
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


		//if (len==1) break;
	}	
	fclose(fions);
	fclose(fin);
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


	printf("totial len=%d\n", len);
	printf("totial Step=%d\n", step);
	printf("Averagge lifetime %lf\n", slift_time/nconst);
	printf("Output ave info\n");
	fprintf(flife, "%lf ", slift_time/nconst);
	fprintf(flife, "%lf\n", avenc);
	fclose(flife);


	return 0;
}



