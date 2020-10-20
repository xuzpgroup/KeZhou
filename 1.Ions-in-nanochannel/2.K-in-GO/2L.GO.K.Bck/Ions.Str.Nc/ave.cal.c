#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_L 20000
#define SAMPLES 4
#define STD 0


int main()
{
	int i,j,st;

	char 	header_char[300];
	char 	filename[300];
	double	ti, data[SAMPLES];
	int 	flag[SAMPLES];
	double 	ss, sq;
	int 	totalsamples, len;
	FILE 	*fin;
	FILE 	*fout;
	
	printf("Start .....\n");

	len=0;
	fin=fopen("./His.nc.all.data", "r");
	fout=fopen("His.nc.ave.dat", "w");
	//printf("%ld\n", len);
	// fgets(header_char, 300, fin);
	// puts(header_char);
	
	for (i=0;i<SAMPLES;i++) flag[i]=1;

	// selected samples 
	// if don't want to use the data, just set flag to 0
	//flag[3]=0;
	//flag[7]=0;
	totalsamples=0;
	for (i=0;i<SAMPLES;i++)
	{
		totalsamples+=flag[i];
	}
	ss=0.0;
	sq=0.0;
	//fprintf(fout, "%lf  %lf  %lf  %lf  %lf\n", 0.0, 0.0, 0.0, 0.0, 0.0);
	fgets(header_char, 300, fin);
	fgets(header_char, 300, fin);
	while (feof(fin)==0)
	{
		len++;
		ss=0.0;
		fscanf(fin, "%lf", &ti);
		for (i=0;i<SAMPLES;i++)
		{
			fscanf(fin, "%lf", &data[i]);
			ss+=data[i]*flag[i];			
		}
		// fscanf(fin, "\n");
		ss=ss/totalsamples;
		sq=0.0;
		for (i=0;i<SAMPLES;i++)
		{
			sq+=(data[i]-ss)*(data[i]-ss)*flag[i];			
		}
		sq=sq/totalsamples;
		sq=sqrt(sq);

		fprintf(fout, "%.0lf  %lf", ti, ss);
		if (STD) fprintf(fout, "%lf  %lf", ss-sq, ss+sq);
		fprintf(fout, "\n");
	}	



	fclose(fout);
	printf("totial len=%d\n", len);
	fclose(fin);
	return 0;
}



