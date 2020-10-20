#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG_1
#define N_MOL 432
#define MAX_T 20001


int main()
{
	int len;
	int bins, step;
	int i,j,st;

	char header_char[300];
	char filename[300];
	double data[8];

	FILE *fin;
	FILE *fout;
	printf("Start .....\n");

	len=0;
	//dumpxu=fopen("mol.dump.xu.xyz", "r");
	fin=fopen("./rdf.acm.data", "r");
	//fout=fopen("rdf.reduced.dat", "w");
	//printf("%ld\n", len);
	fgets(header_char, 300, fin);
	puts(header_char);
	fgets(header_char, 300, fin);
	fgets(header_char, 300, fin);
	

	while (feof(fin)==0)
	{
		len++;
		fscanf(fin, "%d%d\n", &step, &bins);
		sprintf(filename, "./RDF/rdf.%d.data", step);
		puts(filename);
		fout=fopen(filename, "w");
		fputs(header_char, fout);
		// printf("%d %d\n", step, bins);
		for (i=0;i<bins;i++)
		{
			for (j=0;j<6;j++) fscanf(fin, "%lf", &data[j]);
			fscanf(fin, "\n");
			for (j=0;j<6;j++) fprintf(fout, "%12.6lf ", data[j]);
			fprintf(fout, "\n");			
		}
		fprintf(fout, "\n");
		fclose(fout);
	}	

	printf("totial len=%d\n", len);
	fclose(fin);
	return 0;
}



