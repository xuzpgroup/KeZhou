#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define MAXNATOM 20000

int		id[MAXNATOM];
int 	itype[MAXNATOM];
double 	r[MAXNATOM][3];
char 	ele[28][10];
double 	mass[28];

void swap(int i, int j)
{
    int temp;
    double tempdou;
    int kk;

    temp=id[i];		id[i]=id[j];		id[j]=temp;
	temp=itype[i];	itype[i]=itype[j];	itype[j]=temp;
	for (kk=0; kk<3; kk++)
	{
		tempdou=r[i][kk];	
		r[i][kk]=r[j][kk];
		r[j][kk]=tempdou;
	}
}


void quicksort(int left, int right)
{
    int pivot, i, j;

    if (left >= right) { return; }
    pivot = id[left];
    i = left + 1;
    j = right;
    while (1)
    {
        while (i <= right)
        {
            if (id[i] > pivot)
            {
                break;
            }
            i++;
        }
        while (j > left)
        {
            if (id[j] < pivot)
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


void Modify_Element()
{
	sprintf(ele[0], "C");	
	sprintf(ele[1], "C");	
	sprintf(ele[2], "O");	
	sprintf(ele[3], "C");
	sprintf(ele[4], "O");
	sprintf(ele[5], "H");
	sprintf(ele[6], "H");
	sprintf(ele[7], "H");
	sprintf(ele[8], "H");
	sprintf(ele[9], "H");
	sprintf(ele[10], "C");
	sprintf(ele[11], "C");
	sprintf(ele[12], "C");
	sprintf(ele[13], "N");
	sprintf(ele[14], "Cl");
	sprintf(ele[15], "F");
	sprintf(ele[16], "Al");
	sprintf(ele[17], "B");
	sprintf(ele[18], "O");
	sprintf(ele[19], "C");

}
void Modify_Mass()
{
	mass[0]=12.0;		// CA
	mass[1]=12.0;		// CT
	mass[2]=16.0;		// OS
	mass[3]=12.0;		// CF
	mass[4]=16.0;		// OH
	mass[5]=1.0;		// HO

	mass[6]=1.0;		// HC
	mass[7]=1.0;		// H1
	mass[8]=1.0;		// H4
	mass[9]=1.0;	 	// H5
	mass[10]=12.0;		// CT
	mass[11]=12.0;		// CR -
	mass[12]=12.0;		// CW -
	mass[13]=14.0;		// NA -
	mass[14]=35.45;		// Cl
	mass[15]=19.0;		// F
	mass[16]=27.0;		// Al
	mass[17]=10.81;		// B
	mass[18]=10.81;		// O-co2
	mass[19]=10.81;		// C-co2


}

int main()
{
	int len;
	int step, natoms, nmol;
	int i, j, st, tp;
	double selected;
	double totm, totm5;

	char header_char[300];
	char filename[300];
	double data[8];
	double mc[3], m5[3];


	FILE *fin;
	FILE *fout3d;
	FILE *fout2d;
	FILE *fmc;
	FILE *fm5;

	printf("Start .....\n");
	Modify_Element();
	Modify_Mass();

	selected=0.2;

	len=0;

	fin=fopen("../IL.unwrap.xyz", "r");
	fout3d=fopen("../IL.unwrap.3d.xyz", "w");	
	fout2d=fopen("../IL.unwrap.2d.xyz", "w");
	fmc=fopen("../IL.mc.xyz", "w");
	fm5=fopen("../IL.m5.xyz", "w");
	
	while (feof(fin)==0)
	{
		len++;
		//if (len % 10==0) printf("Now step = %d", step);
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &step);
		fgets(header_char, 300, fin);
		fscanf(fin, "%d\n", &natoms);
		for (i=0;i<5;i++) fgets(header_char, 300, fin);
	
		if (len % 100 ==0) printf("Now is %d\n", step);
		fprintf(fout2d, "%d\n", natoms);
		fprintf(fout3d, "%d\n", natoms);
		fprintf(fout2d, "Atoms. Timestep: %d\n", step);
		fprintf(fout3d, "Atoms. Timestep: %d\n", step);
		for (i=0;i<natoms;i++)
		{
			fscanf(fin, "%d%d", &id[i], &itype[i]);
			fscanf(fin, "%lf%lf%lf\n", &r[i][0], &r[i][1], &r[i][2]);
		}
		//quicksort(0, natoms-1);
		for (i=0;i<natoms;i++)
		{
			fprintf(fout3d, "%s  ", ele[itype[i]-1+6]);
			// fprintf(fout3d, "%d  ", id[i]);
			fprintf(fout3d, "%12.6lf%12.6lf%12.6lf\n", r[i][0], r[i][1], r[i][2]);	
			fprintf(fout2d, "%s  ", ele[itype[i]-1]);
			fprintf(fout2d, "%12.6lf%12.6lf%12.6lf\n", r[i][0], r[i][1], 0.0);
		}

		nmol=natoms/30;

		fprintf(fmc, "%d\n", nmol*2);
		fprintf(fmc, "Atoms. Timestep: %d\n", step);
		fprintf(fm5, "%d\n", nmol*2);
		fprintf(fm5, "Atoms. Timestep: %d\n", step);
	

		// for BMIM
		for (i=0;i<nmol;i++)
		{
			st=i*25;
			mc[0]=0.0;			m5[0]=0.0;
			mc[1]=0.0;			m5[1]=0.0;
			mc[2]=0.0;			m5[2]=0.0;

			totm=0.0;
			totm5=0.0;

			for (j=st; j<st+25; j++)
			{
				tp=itype[j]-1+6;
				mc[0]+=r[j][0]*mass[tp];
				mc[1]+=r[j][1]*mass[tp];
				mc[2]+=r[j][2]*mass[tp];
				totm=totm+mass[tp];
				// if (itype[j] ==12 || itype[j] ==13 || itype[j] ==14)
				if (tp+1 ==12 || tp+1==13 || tp+1==14)
				{
					m5[0]+=r[j][0]*mass[tp];
					m5[1]+=r[j][1]*mass[tp];
					m5[2]+=r[j][2]*mass[tp];
					totm5=totm5+mass[tp];					
				}
			}		
			// printf("%d %lf\n", i+1, totm);
			mc[0]=mc[0]/totm;
			mc[1]=mc[1]/totm;
			mc[2]=mc[2]/totm;
			m5[0]=m5[0]/totm5;
			m5[1]=m5[1]/totm5;
			m5[2]=m5[2]/totm5;
			
			fprintf(fmc, "C  ");
			fprintf(fmc, "%12.6lf%12.6lf%12.6lf\n", mc[0], mc[1], mc[2]);
			// fprintf(fmc, "%12.6lf%12.6lf%12.6lf\n", mc[0], mc[1], 0.0);
			fprintf(fm5, "C  ");
			fprintf(fm5, "%12.6lf%12.6lf%12.6lf\n", m5[0], m5[1], m5[2]);
			// fprintf(fm5, "%12.6lf%12.6lf%12.6lf\n", m5[0], m5[1], 0.0);

		}

		// for BF4
		for (i=0;i<nmol;i++)
		{
			st=nmol*25+i*5;
			mc[0]=0.0;
			mc[1]=0.0;
			mc[2]=0.0;
			totm=0.0;
			for (j=st; j<st+5; j++)
			{
				tp=itype[j]-1+6;
				mc[0]+=r[j][0]*mass[tp];
				mc[1]+=r[j][1]*mass[tp];
				mc[2]+=r[j][2]*mass[tp];
				totm+=mass[tp];
			}
			mc[0]/=totm;
			mc[1]/=totm;
			mc[2]/=totm;
			// printf("%d %lf\n", i+1, totm);
			fprintf(fmc, "B  ");
			fprintf(fmc, "%12.6lf%12.6lf%12.6lf\n", mc[0], mc[1], mc[2]);
			fprintf(fm5, "B  ");
			fprintf(fm5, "%12.6lf%12.6lf%12.6lf\n", mc[0], mc[1], mc[2]);
		}
		// for gas 
		// fprintf(fmc, "N  ");
		// fprintf(fmc, "%12.6lf%12.6lf%12.6lf\n", r[natoms-1][0], r[natoms-1][1], r[natoms-1][2]);
		// fprintf(fm5, "N  ");
		// fprintf(fm5, "%12.6lf%12.6lf%12.6lf\n", r[natoms-1][0], r[natoms-1][1], r[natoms-1][2]);

		//if (len==1) break;
	}	

	fclose(fout2d);
	fclose(fout3d);
	fclose(fmc);
	fclose(fm5);
	
	printf("totial len=%d\n", len);
	printf("totial Step=%d\n", step);

	return 0;
}



