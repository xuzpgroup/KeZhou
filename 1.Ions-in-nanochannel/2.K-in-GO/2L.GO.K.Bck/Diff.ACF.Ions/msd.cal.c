#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_L 20000
#define SAMPLES 5


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
	fin=fopen("./msd.all.data", "r");
	fout=fopen("msd.ave.dat", "w");
	fgets(header_char, 300, fin);
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
	fprintf(fout, "%lf  %lf  %lf  %lf  %lf\n", 0.0, 0.0, 0.0, 0.0, 0.0);
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

		// diffusion unit 10^-9m2/s
		if (len % 10==0) fprintf(fout, "%lf  %lf  %lf  %lf  %lf\n", ti, ss, ss/ti/4.0*10, ss-sq, ss+sq);
	}	



	fclose(fout);
	printf("totial len=%d\n", len);
	fclose(fin);
	return 0;
}



