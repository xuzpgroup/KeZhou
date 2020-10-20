#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "fftw3.h"

#define CAL_MACF
#define MAXNATOM 120
#define MAX_CLU 30
//#define MAXT 131072
#define MAXT 65536

// #define MAXT 1048576	 // 10ns each 10 fs
// #define MAXT  524288  // 5ns each 10 fs

#define OUTSCALE 5
#define REDUCED 2

// ############## for bilayer
#define R_TS 3.58400
#define R_CIP 2.84800
#define R_SIP 4.91200
#define R_TS2 6.0


short MM[MAXNATOM][MAXNATOM][MAXT];

double 	rM[MAXNATOM][3];
double 	rC[MAXNATOM][3];
double  dist[MAXNATOM][MAXNATOM];
double  autof[MAXT];
double	Macf[MAXT];

int len;
int step, natoms;


/*Union-Find*/
/*##########################*/
int pre[MAXNATOM*2];
int top[MAX_CLU][MAX_CLU];
int ptp[MAXNATOM*2][2];     // the number of ion in cluster i
int flag_free[MAXNATOM*2];  // is the ion i is free ion

double stop[MAX_CLU][MAX_CLU];

/*##########################*/


// for histogram
// int    bins;
// double hst[HisItems], hed[HisItems], delt[HisItems];
// double his[MAXBINS][HisItems];
// 	// 	info
// 	// 	0	theta 1
// 	// 	1	cos (theta1)


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

int log_2(int s)
{
	int i, ss;
	i=1;
	ss=s;
	while (ss>1)
	{
		ss=ss/2;	i++;
	}
	return i;
}

int exp_2(int s)
{
	int i, ss;
	i=s;
	ss=1;
	while (i>0)
	{
		ss=ss*2;	i--;
	}
	return ss;
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



double Vec_Dot(double *a, double *b)
{
	double s;
	s=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	return s;
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

 
int unionsearch(int root)
{
	int son, tmp;
	son = root;
	while(root != pre[root]) //寻找掌门ing……
		root = pre[root];
	while(son != root) //路径压缩
	{
		tmp = pre[son];
		pre[son] = root;
		son = tmp;
	}
	return root; //掌门驾到~
}
 
void clear_top(int nat)
{
	int i,j;
	// construct the toplogy
	for (i=0;i<MAX_CLU;i++)
	{
		for (j=0;j<MAX_CLU;j++)
		{
			top[i][j]=0;
		}
	}

	for (i=0;i<nat*2;i++)
	{
		ptp[i][0]=0;
		ptp[i][1]=0;	
		flag_free[i]=0; // is free ions or not		
	}

}

int Top_Ionic_Net(int ti) // return number of free ions
{
	int i, j, start, end, root1, root2;
	int clut;
	int root, nan, nct;
	int nfree;

	for (i=0;i<natoms*2;i++) pre[i]=i;
	for (i=0;i<natoms;i++)
	{
		for (j=0;j<natoms;j++)
		{
			if (MM[i][j][ti]==1)
			{
				root1 = unionsearch(i);
				root2 = unionsearch(j+natoms);
				if(root1 != root2) //掌门不同？踢馆！~
				{
					pre[root1] = root2;
				}				
			}
		}
	}
	clear_top(natoms);

	for (i=0;i<natoms*2;i++)
	{
		root=unionsearch(i);
		if (i<natoms)
		{
			ptp[root][0]++;
		}
		else 
		{
			ptp[root][1]++;
		}
	}
	for (clut=0;clut<natoms*2;clut++)
	{
		nct=ptp[clut][0];
		nan=ptp[clut][1];
		if (nan+nct>0.55)
		{
			top[nct][nan]++;
		}		
	}

	nfree=0; // number of free ions
	for (i=0;i<natoms*2;i++)
	{
		root=unionsearch(i);
		if (ptp[root][0]+ptp[root][1]==1)
		{
			flag_free[i]=1;
			nfree++;
		}
	}

	for (i=0;i<MAX_CLU;i++)
	{
		for (j=0;j<MAX_CLU;j++)
		{
			stop[i][j]+=top[i][j];
		}
	}

	return nfree;
}



// void HistogramIt(double samples, int jj)
// {
// 	int kk; 

// 	kk=(samples-hst[jj])/delt[jj];
// 	if (kk<bins)
// 	{
// 		his[kk][jj]+=1.0;
// 	}
// }

int main()
{
	int i, j, st, tp, k;
	int id, itype, imol, tt;
	double chg, vx, vy, vz, xu, yu, zu;
    double xlo, xhi, ylo, yhi, zlo, zhi;
	double box[3];
	double timestep;
	double v;
	char header_char[300];
	char filename[300];
	int free_ion;

	double density_pair_cip, density_pair_sip;
	double ave_density_pair_cip, ave_density_pair_sip;

	FILE *fin;
	FILE *fpt;
	FILE *fcluster;
	FILE *ftop;
	
//  for FFTW
    fftw_complex *in, *out;
    fftw_plan p;

	box[0]=32.078416;
	box[1]=32.078416;
	box[2]=32.078416;

	// important paramters
	timestep=0.010*REDUCED; // ps
	
	printf("Start .....\n");

	len=0;
	tt=0;

	ave_density_pair_cip=0.0;
	ave_density_pair_sip=0.0;

	fin=fopen("../ions.dat", "r");
	fpt=fopen("./pair.t.dat", "w");
	fcluster=fopen("./cluster.dat", "w");
	
	
	while (feof(fin)==0)
	{
		tt++;
		// if (tt>800000) break;
		if (tt % REDUCED == 1) len++;
		// len++;
		if (len % 1000 ==1) printf("Now is %d\n", len);
		if (len>MAXT) printf("ERROR!!!\n");
		// printf("time = %d\n", tt);	
		// input water and ions from FILE: ions.dat
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &step);
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &natoms);
		for (i=0;i<1;i++) fgets(header_char, 300, fin);
	    fscanf(fin, "%lf%lf\n", &xlo, &xhi);
	    fscanf(fin, "%lf%lf\n", &ylo, &yhi);
	    fscanf(fin, "%lf%lf\n", &zlo, &zhi);
	    fgets(header_char, 300, fin); 
	    
	    box[0]=xhi-xlo;
	    box[1]=yhi-ylo;
	    box[2]=zhi-zlo;
	    
		natoms=natoms/2;
		// printf(" Natom = %d \n", natoms);
		for (i=0;i<natoms*2;i++)
		{
			// ITEM: ATOMS id type xu yu zu
			fscanf(fin, "%d%d%d%lf", &id, &imol, &itype, &chg);
			fscanf(fin, "%lf%lf%lf", &vx, &vy, &vz);
			fscanf(fin, "%lf%lf%lf\n", &xu, &yu, &zu);
			if (itype==9)
			{
				rM[i][0]=xu;
				rM[i][1]=yu;
				rM[i][2]=zu;
			}	
			if (itype==10)
			{
				rC[i-natoms][0]=xu;
				rC[i-natoms][1]=yu;
				rC[i-natoms][2]=zu;
			}					
		}
		if (tt % REDUCED != 1) continue;

		// calculate distance between ions
		// get MM matrix
		// printf("calculate distance between ions\n");
		for (i=0;i<natoms;i++)
		{
			for (j=0;j<natoms;j++)
			{
				dist[i][j]=Distant_PBC(rM[i], rC[j], box);
				MM[i][j][len-1]=0;
				if (dist[i][j]<R_TS)
				{
					MM[i][j][len-1]=1;
				}
				else if (dist[i][j]<R_TS2)
				{
					MM[i][j][len-1]=2;
				}
			}
		}

		// calculate pair density 
		density_pair_cip=0.0;
		density_pair_sip=0.0;
		for (i=0;i<natoms;i++)
		{
			for (j=0;j<natoms;j++)
			{
				if (MM[i][j][len-1]==1)
				{
					density_pair_cip=density_pair_cip+1.0;
				}
				if (MM[i][j][len-1]==2)
				{
					density_pair_sip=density_pair_sip+1.0;
				}	
			}
		}
		density_pair_cip=density_pair_cip/natoms;
		density_pair_sip=density_pair_sip/natoms;

		fprintf(fpt, "%12.6lf  %lf  %lf  %lf  %lf \n", \
			timestep*(len-1), \
			density_pair_cip, density_pair_cip*natoms, \
			density_pair_sip, density_pair_sip*natoms);
		
		ave_density_pair_sip+=density_pair_sip;
		ave_density_pair_cip+=density_pair_cip;

		free_ion=Top_Ionic_Net(len-1);
		// output cluster
		// printf("Num_culster = %d\n", natoms*2-free_ion);

		fprintf(fcluster, "ITEM: TIMESTEP\n");
		fprintf(fcluster, "%d\n", len);
		fprintf(fcluster, "ITEM: NUMBER OF ATOMS\n");
		fprintf(fcluster, "%d\n", natoms*2-free_ion);
		fprintf(fcluster, "ITEM: BOX BOUNDS pp pp pp\n");
		fprintf(fcluster, "%lf %lf\n", 0.0, box[0]);
		fprintf(fcluster, "%lf %lf\n", 0.0, box[1]);
		fprintf(fcluster, "%lf %lf\n", 0.0, box[2]);
		fprintf(fcluster, "ITEM: ATOMS id type q xu yu zu \n");

		j=0;
		for (i=0;i<2*natoms;i++)
		{
			if (i<natoms) 
			{
				itype=9;
				chg=1.0;
				xu=rM[i][0];
				yu=rM[i][1];
				zu=rM[i][2];	
			}
			else 
			{
				itype=10;
				chg=-1.0;
				xu=rC[i-natoms][0];
				yu=rC[i-natoms][1];
				zu=rC[i-natoms][2];					
			}
			if (!flag_free[i])
			{
				j++;
				fprintf(fcluster, "%d ", j);
				fprintf(fcluster, "%d ", itype);
				fprintf(fcluster, "%lf ", chg);
				fprintf(fcluster, "%lf %lf %lf \n", xu, yu, zu);
			}
		}
	
		// if (len>=1) break;
	}	
	
	printf(" The total lenngth of data is %d \n", len);
	
	double check_top;
	check_top=0.0;
	ftop=fopen("./cluster.top.dat", "w");
	
	for (i=0;i<MAX_CLU;i++)
	{
		for (j=0;j<MAX_CLU;j++)
		{
			stop[i][j]/=len;
			fprintf(ftop, "%d %d %lf\n", i, j, stop[i][j]);
			check_top+=stop[i][j]*(i+j);
		}
		fprintf(ftop, "\n");
	}
	printf("check_top = %lf\n", check_top);
	fclose(ftop);


	ave_density_pair_sip/=len;
	ave_density_pair_cip/=len;
	fprintf(fpt, "####  %lf  %lf  %lf  %lf \n", \
		ave_density_pair_cip, ave_density_pair_cip*natoms,\
		ave_density_pair_sip, ave_density_pair_sip*natoms);
	
	// fclose(fions);
	fclose(fin);
	fclose(fpt);

#ifdef CAL_MACF
	int total_l;
    FILE *f0, *f2;

	printf(" Calculate M acf by FFT ...\n");
    printf("Total step = %d\n", tt);
    printf("Selected step = %d\n", len);

    // total_l=exp_2(log_2(len*2-1));
    // total_l=exp_2(log_2(len/2-1));
    total_l=MAXT;

    printf("FFT Total length = %d\n", total_l);

    for (i=0;i<natoms;i++)
    {
    	for (j=0;j<natoms;j++)
    	{
	        for (k=len;k<total_l;k++)
	        {
	            MM[i][j][k]=0;
	            MM[i][j][k]=0;
	            MM[i][j][k]=0;
	        }    		
    	}
      
    }
    for (i=0;i<total_l;i++)
    {
        Macf[i]=0.0;
    }
    printf("\n Start .... \n");
	for (i=0;i<natoms;i++)
	{
		printf(" >>>>> FFT for %d x j \n", i+1);
		for (j=0;j<natoms;j++)
		{
	        // FFT
	        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total_l);
	        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total_l);
	        //printf("Forward FFT %d \n", k);
	        for (k=0;k<total_l;k++)
	        {
	            // in[k][0]=(MM[i][j][k]==1);
	            in[k][0]=0.0;
	            if (MM[i][j][k]==1) in[k][0]=1.0;
	            in[k][1]=0.0;
	        }
	        p = fftw_plan_dft_1d(total_l, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	        fftw_execute(p);
	        fftw_destroy_plan(p);
	            
	        // Inverse FFT
	        // printf("Forward FFT\n");
	        for (k=0;k<total_l;k++)
	        {
	            v=out[k][0]*out[k][0]+out[k][1]*out[k][1];
	            in[k][0]=v;
	            in[k][1]=0.0;
	        }

	        //printf("Backward FFT %d \n", k);
	        p = fftw_plan_dft_1d(total_l, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	        fftw_execute(p);
	        fftw_destroy_plan(p);

	        for (k=0;k<total_l;k++)
	        {
	            autof[k]=out[k][0]/total_l/(len-k);
	        }
	        fftw_free(in); 
	        fftw_free(out);
	        for (k=0;k<total_l;k++)
	        {
	            Macf[k]+=autof[k];
	        }
		}	
	}

	for (i=0;i<total_l;i++) Macf[i]/=natoms;

    // output ACF 
    printf("Output ACF\n");
    double scaleacf;

    scaleacf=1.0; // unit m/s * m/2
    f2=fopen("./0rtcf.data", "w");
    fprintf(f2, "#### t[ps] ### MACF <V(0)*V(t)> \n");
    for (i=0;i<len/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+0.5)*timestep);
        fprintf(f2, "  %20.6E  ", Macf[i]*scaleacf);
        fprintf(f2, "\n");
    }
    fclose(f2);

    // output normalized ACF 
    f2=fopen("./01_norm_rtcf.data", "w");
    fprintf(f2, "#### t[ps] ### MACF <Q(0)*Q(t)> \n");
    for (i=0;i<len/OUTSCALE;i++)
    {
        fprintf(f2, "%-12.6lf", (i+0.5)*timestep);
        fprintf(f2, "  %20.6E  ", Macf[i]/Macf[0]);
        fprintf(f2, "\n");
    }
    fclose(f2);

#endif


	return 0;
}



