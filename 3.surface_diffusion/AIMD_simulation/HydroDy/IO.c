// ################
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

/*code for input and output*/
#define MAX_ATOMS 500
#define PI 3.1415926


void write_hander()
{
	printf("// #####################################################################\n");
	printf("// #                                                                   #\n");	
	printf("// # This program is used for analysing the hydration shell under nano #\n");	
	printf("// # confinement, including the statical structure and dynamical       #\n");	
	printf("// # infromation.                                                      #\n");	
	printf("// #                                                                   #\n");	
	printf("// #####################################################################\n");	
	printf("\n");
}


void write_tail()
{
	printf("// #####################################################################\n");
	printf("// #                                                                   #\n");	
	printf("// #                      End                                          #\n");	
	printf("// #                                                                   #\n");	
	printf("// #####################################################################\n");	

}

void dump_into_box(double *box, char *filename)
{
    int    DIM;
    double r[MAX_ATOMS][3];
    int    tp[MAX_ATOMS];
    char   char_tp;
    int    i, len, natoms, total_len;
    int    discard_frames;
 
    char header_char[300];
    char output_file[300];

    FILE *fin;
    FILE *fout;

    printf(" Dump into box \n");

    len=0;   
    total_len=0;
    discard_frames=0;


    fin=fopen(filename, "r");
	
	sprintf(output_file, "dump.dat");
    fout=fopen(output_file, "w");  
    printf(" >> Output wrap coordinate, filename: %s\n", output_file);
    
    while (feof(fin)==0)
    {
        len++;
        if (len % 1000==0) printf("Now is %d\n", len);
        fscanf(fin, "%d\n", &natoms);
        fgets(header_char, 300, fin);
        // puts(header_char);
        // printf("Natoms = %d \n", natoms);
        for (i=0;i<natoms;i++)
        {
            fscanf(fin, "  %c%lf%lf%lf\n", &char_tp, \
                &r[i][0], &r[i][1], &r[i][2]);
            if (char_tp=='C') tp[i]=1;
            if (char_tp=='H') tp[i]=7;
            if (char_tp=='O') tp[i]=8;
        }
      
		fprintf(fout, "ITEM: TIMESTEP\n");       
		fprintf(fout, "%d\n", len);       
		fprintf(fout, "ITEM: NUMBER OF ATOMS\n");       
		fprintf(fout, "%d\n", natoms);
		fprintf(fout, "ITEM: BOX BOUNDS pp pp pp\n");    
		fprintf(fout, "%lf %lf \n", 0.00, box[0]);
		fprintf(fout, "%lf %lf \n", 0.00, box[1]);
		fprintf(fout, "%lf %lf \n", 0.00, box[2]);
		fprintf(fout, "ITEM: ATOMS id type xu yu zu \n");    
		for (i=0;i<natoms;i++)
		{
			fprintf(fout, "%d ", i+1);
			fprintf(fout, "%d ", tp[i]);
			fprintf(fout, "%lf ", r[i][0]);
			fprintf(fout, "%lf ", r[i][1]);
			fprintf(fout, "%lf ", r[i][2]);
			fprintf(fout, "\n");
		}
        total_len++;   
        // printf(" check %d, %lf\n", len, his[40]);  
    }

}