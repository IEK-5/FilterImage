#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include "floatimage.h"
#include "filter.h"

typedef enum {DX, DY, NABLA, PEAK} filter_type;

void Help(struct option OPTS[], char *desc[])
{
	int i=0;
	printf("Usage: FilterImage [options] -i input-image [-o/-t] output-image\n");
	printf("Options:\n");
	while (OPTS[i].name)
	{
		if (isalnum(OPTS[i].val))
			printf("--%s [-%c]", OPTS[i].name, OPTS[i].val);
		else
			printf("--%s     ", OPTS[i].name);
		if (desc[i])
			printf(":\t%s\n", desc[i]);
		else
			printf("\n");
		
		i++;
	}
}

//#define DEBUG	
int main(int argc, char **argv)
{
	image Iin, Iout;
	char *fin=NULL, *fout=NULL, *ftxt=NULL;
	int N=3, m=3; 
	filter_type T=NABLA;
	int deriv_order=0, ItoV=0, FFT=0, p=1, c;
	double Temp; 
	char d='n';
	
	
	
	while (1)
	{
		static struct option long_options[] =
		{
			{"in",          required_argument, 0, 'i'},
			{"out",         required_argument, 0, 'o'},
			{"txt-out",     required_argument, 0, 't'},
			{"I2V",         required_argument, 0, 'V'},
			{"N-filter",    required_argument, 0, 'N'},
			{"filter-order",required_argument, 0, 'm'},
			{"dx",          required_argument, 0, 'x'},
			{"dy",          required_argument, 0, 'y'},
			{"nabla",       required_argument, 0, 'n'},
			{"peak",        required_argument, 0, 'p'},
			{"fft",               no_argument, 0, 'f'},
			{"help",              no_argument, 0, 'h'},
			{0, 0, 0, 0}
		};
		char *desc[]={
			"\tspecify input-image file name (arg=filename)",
			"\tspecify output-image filename (arg=filename)",
			"\tspecify output text image filename\n\t\t\t(arg=filename or - for stdout)",
			"\tconvert luminescence intensity to junction voltage\n\t\t\t(arg=temperature in K)",
			"Size of the filter (arg=integer>0)",
			"filter-order (arg=integer>0)",
			"\tcreate derivative to x (arg=integer,\n\t\t\t0<=arg<=filter-order, the order of the derivative)",
			"\tcreate derivative to y (arg=integer,\n\t\t\t0<=arg<=filter-order, the order of the derivative)",
			"\tnabla operator (dx+dy) (arg=integer,\n\t\t\t0<=arg<=filter-order, the order of the derivative)",
			"\tpeak detector (arg=integer>0 peak size)",
			"\tuse the Fast Fourier Transform, makes computation\n\t\t\ttime independent of filter size but may lead to\n\t\t\tedge effects at the borders of the image",
			"\tshow this help message and exit",
			NULL};
			
		/* getopt_long stores the option index here. */
		int option_index = 0;
		c = getopt_long (argc, argv, "i:o:V:N:m:x:y:n:p:ft:h",long_options, &option_index);
		/* Detect the end of the options. */
		if (c == -1)
			break;
		switch (c)
		{
			case 'h':
				Help(long_options,desc);
				exit(0);
			case 'i':
				if (!optarg)
				{
					fprintf(stderr, "Error: --in requires an input file\n");
					return 1;	
				}
				
				fin=malloc((strlen(optarg)+1)*sizeof(char));
				strcpy(fin, optarg);
				break;
			case 'o':
				if (!optarg)
				{
					fprintf(stderr, "Error: --out requires an output file\n");
					return 1;	
				}
				fout=malloc((strlen(optarg)+1)*sizeof(char));
				strcpy(fout, optarg);
				break;
			case 't':
				if (!optarg)
				{
					fprintf(stderr, "Error: --out requires an output file\n");
					return 1;	
				}
				ftxt=malloc((strlen(optarg)+1)*sizeof(char));
				strcpy(ftxt, optarg);
				break;
			case 'V':
				if (!optarg)
				{
					fprintf(stderr, "Error: --out requires an output file\n");
					return 1;	
				}
				ItoV=1;
				Temp=atof(optarg);
				break;
			case 'N':
				if (!optarg)
				{
					fprintf(stderr, "Error: --L_filter requires a length specification for the filter\n");
					return 1;	
				}
				N=atoi(optarg);
				break;
			case 'm':
				if (!optarg)
				{
					fprintf(stderr, "Error: --L_filter requires a length specification for the filter\n");
					return 1;	
				}
				m=atoi(optarg);
				break;
			case 'x':
				if (!optarg)
				{
					fprintf(stderr, "Error: --Ident requires an epsilon specification for the filter (default 1e-12)\n");
					return 1;	
				}
				T=DX;
				deriv_order=atoi(optarg);
				d='x';
				break;
			case 'y':
				if (!optarg)
				{
					fprintf(stderr, "Error: --Ident requires an epsilon specification for the filter (default 1e-12)\n");
					return 1;	
				}
				T=DY;
				deriv_order=atoi(optarg);
				d='y';
				break;
			case 'n':
				if (!optarg)
				{
					fprintf(stderr, "Error: --Ident requires an epsilon specification for the filter (default 1e-12)\n");
					return 1;	
				}
				T=NABLA;
				deriv_order=atoi(optarg);
				d='n';
				break;
			case 'p':
				if (!optarg)
				{
					fprintf(stderr, "Error: --Ident requires an epsilon specification for the filter (default 1e-12)\n");
					return 1;	
				}
				T=PEAK;
				p=atoi(optarg);
				break;
			case 'f':
				FFT=1;
				break;
			case '?':
				/* getopt_long already printed an error message. */
				break;
			default:
				abort ();
			}
		}
	if (!fin)
	{
		fprintf(stderr,"Error: No input file\n");
		exit(1);
	}
	if (!fout && !ftxt)
	{
		fprintf(stderr,"Error: No output specified, will not do work for nothing\n");
		exit(1);
	}
	if (strncmp(fin, "--test-image--", 15)==0)
		Iin=TestImage(100,100,2);
	else
		Iin=FloatimageRead(fin);
	free(fin);
	
	if (ItoV)
		EL2Vj(Iin, Temp); // convert to junction voltages
	
	switch(T)
	{
		case DX:
		case DY:
		case NABLA:
			if (FFT)
				Iout=FFT_PolynomalFilter(Iin, N, m, deriv_order, d); // apply the filter
			else
				Iout=FFT_PolynomalFilter(Iin, N, m, deriv_order, d); // apply the filter
			break;
		case PEAK:
			if (FFT)
				Iout=FFT_PolynomalExtremaLocator(Iin, N, p);
			else
				Iout=PolynomalExtremaLocator(Iin, N, p);
			break;
		default:
			fprintf(stderr, "Error: check status of nuclear power planst in the vicinity cause this error cannot happen\n");
			exit(1);
	}
	if (fout)
	{
		FloatimageWrite(fout, Iout, 1,0,1);
		free(fout);
	}
	
	if (ftxt)
	{
		if (strncmp(fin, "-", 15)==0)
			Floatimage2stdout(Iout);
		else
			FloatimageTXTWrite(ftxt, Iout);
		free(ftxt);
	}		
		
	FreeImage(&Iin);
	FreeImage(&Iout);
	return 0;
}
