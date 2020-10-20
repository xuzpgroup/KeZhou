#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define Height 6.8


int main()
{
	int ii, step;
	double x,y,z; 
	char header_char[300];
	char filename[300];
	
	FILE *fin, *fout;

	printf("Start .....\n");

	fout=fopen("pos.data", "w");

	for (ii=0;ii<1;ii++)
	{
		sprintf(filename, "../1.H6.8.M1.Na.temp.%d/pos.pmf.data", ii+1);\
		puts(filename);
		fin=fopen(filename, "r");
		fgets(header_char, 300, fin);
		fgets(header_char, 300, fin);
		while (feof(fin)==0)
		{
			fscanf(fin, "%d%lf%lf%lf\n", &step, &x, &y, &z);
			z=z-Height/2.0;
			if (z<0) z+=Height*2;

			if (z>Height*2) z=z-2*Height;
			// if (z>Height) z=2*Height-z;
			
			fprintf(fout, "%d %lf %lf %lf\n", step, x, y, z);
			fprintf(fout, "%d %lf %lf %lf\n", step, x, y, 2*Height-z);	
		}
		fclose(fin);	
	}
	fclose(fout);

	return 0;
}



