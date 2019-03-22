#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include "floatimage.h"
#include "floatimage_io.h"
#include "filter.h"

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
	int Nx=3, Ny=3, m=3, stepx=1, stepy=1; 
	int deriv_order=0, ItoV=0, FFT=1, c;
	double Temp; 
	double fx=1;
	double fy=1;
	
	
	
	while (1)
	{
		static struct option long_options[] =
		{
			{"in",          required_argument, 0, 'i'},
			{"out",         required_argument, 0, 'o'},
			{"txt-out",     required_argument, 0, 't'},
			{"I2V",         required_argument, 0, 'V'},
			{"N-filter",    required_argument, 0, 'N'},
			{"Nxy-filter",  required_argument, 0, 'R'},
			{"step",        required_argument, 0, 's'},
			{"stepxy",      required_argument, 0, 'r'},
			{"polynomal-order",required_argument, 0, 'm'},
			{"filter-order",required_argument, 0, 'M'},
			{"dx",          required_argument, 0, 'x'},
			{"dy",          required_argument, 0, 'y'},
			{"nabla",             no_argument, 0, 'n'},
			{"fft",               no_argument, 0, 'f'},
			{"plain",             no_argument, 0, 'p'},
			{"help",              no_argument, 0, 'h'},
			{0, 0, 0, 0}
		};
		char *desc[]={
			"\tspecify input-image file name (arg=filename)",
			"\tspecify output-image filename (arg=filename)",
			"\tspecify output text image filename\n\t\t\t(arg=filename or - for stdout)",
			"\tconvert luminescence intensity to junction voltage\n\t\t\t(arg=temperature in K)",
			"Size of the filter (arg=integer>0)",
			"Size of the filter in x and y direction\n\t\t\t(arg=integer,integer>0)",
			"\tstep size for sparse filters (arg=integer>0)",
			"\tstep size in x and y direction for sparse filters\n\t\t\t(arg=integer,integer>0)",
			"polynomal-order (arg=integer>0)",
			"filter-order (arg=integer>0 and <= polynomal-order)",
			"\tcontribution in x direction (arg=float)",
			"\tcontribution in y direction (arg=float)",
			"\tnabla operator (dx+dy)",
			"\tuse the Fast Fourier Transform, makes computation\n\t\t\ttime independent of filter size but may lead to\n\t\t\tedge effects at the borders of the image",
			"\tshow this help message and exit",
			NULL};
			
		/* getopt_long stores the option index here. */
		int option_index = 0;
		c = getopt_long (argc, argv, "i:o:V:N:R:m:M:x:y:ns:r:fpt:h",long_options, &option_index);
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
					fprintf(stderr, "Error: --N_filter requires a length specification for the filter\n");
					return 1;	
				}
				Nx=Ny=atoi(optarg);
				break;
			case 'R':
			{
				char *r;
				if (!optarg)
				{
					fprintf(stderr, "Error: --Nxy_filter requires a length specification for the filter\n");
					return 1;	
				}
				r=optarg;
				while ((*r)&&(*r!=','))
					r++;
				if (!*r)
				{
					fprintf(stderr, "Error: Nxy-filter requires and argument in the form <num>,<num> (without spaces!)\n");
					return 1;
				}
				*r='\0';
				r++;
				Nx=atoi(optarg);
				Ny=atoi(r);
				break;
			}
			case 's':
				if (!optarg)
				{
					fprintf(stderr, "Error: --step requires a length specification for the filter step\n");
					return 1;	
				}
				stepx=stepy=atoi(optarg);
				break;
			case 'r':
			{
				char *r;
				if (!optarg)
				{
					fprintf(stderr, "Error: --stepxy requires a specification for the filter steps\n");
					return 1;	
				}
				r=optarg;
				while ((*r)&&(*r!=','))
					r++;
				if (!*r)
				{
					fprintf(stderr, "Error: --stepxy requires and argument in the form <num>,<num> (without spaces!)\n");
					return 1;
				}
				*r='\0';
				r++;
				stepx=atoi(optarg);
				stepy=atoi(r);
				break;
			}
			case 'm':
				if (!optarg)
				{
					fprintf(stderr, "Error: --polynomal-order requires an order specification\n");
					return 1;	
				}
				m=atoi(optarg);
				break;
			case 'M':
				if (!optarg)
				{
					fprintf(stderr, "Error: --filter-order requires an order specification\n");
					return 1;	
				}
				deriv_order=atoi(optarg);
				break;
			case 'x':
				if (!optarg)
				{
					fprintf(stderr, "Error: --dx requires a magnitude\n");
					return 1;	
				}
				fx=atof(optarg);
				break;
			case 'y':
				if (!optarg)
				{
					fprintf(stderr, "Error: --dy requires a magnitude\n");
					return 1;	
				}
				fy=atof(optarg);
				break;
			case 'n':
				fx=1;
				fy=1;
				break;
			case 'f':
				FFT=2;
				break;
			case 'p':
				FFT=0;
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
		
	switch(FFT)
	{
		case 0:
			Iout=PL_PolynomalFilter(Iin, Ny, Nx, stepy, stepx, m, deriv_order, fx, fy); // apply the filter
			break;
		case 1:
			Iout=PolynomalFilter(Iin, Ny, Nx, stepy, stepx, m, deriv_order, fx, fy); // apply the filter
			break;
		case 2:
			Iout=FFT_PolynomalFilter(Iin, Ny, Nx, stepy, stepx, m, deriv_order, fx, fy); // apply the filter
			break;
		default:
			fprintf(stderr,"Error: the sky is falling and I want my mommy\n");
			exit(1);
	}
	
	if (fout)
	{
		FloatimageWrite(fout, Iout, 1,0,1);
		free(fout);
	}
	
	if (ftxt)
	{
		if (strncmp(fin, "-", 2)==0)
			Floatimage2stdout(Iout);
		else
			FloatimageTXTWrite(ftxt, Iout);
		free(ftxt);
	}		
		
	FreeImage(&Iin);
	FreeImage(&Iout);
	return 0;
}
