#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <wand/MagickWand.h>
#include "error.h"
#include "floatimage.h"

/* magickwand error handler */
void ThrowWandException(MagickWand *wand)
{
	char *description;
	ExceptionType severity;

	description=MagickGetException(wand,&severity);
	(void) fprintf(stderr,"%s %s %lu %s\n",GetMagickModule(),description);
	description=(char *) MagickRelinquishMemory(description);
	exit(-1);
}

image FloatimageRead(char *fn)
{
	image I;
	
	MagickWand *wand = NULL;
	MagickBooleanType status;
    size_t width;
    size_t height;
	
	// init lib
	MagickWandGenesis();
	wand=NewMagickWand();
	
	// open image file
	status=MagickReadImage(wand,fn);
	if (status == MagickFalse)
		ThrowWandException(wand);
		
	status=MagickTransformImageColorspace(wand, GRAYColorspace);
	if (status == MagickFalse)
		ThrowWandException(wand);
	
	
    // get the size and allocate the image
	width = MagickGetImageWidth(wand);
    height = MagickGetImageHeight(wand);
    fprintf(stderr, "image size %lux%lu\n", width, height);
    I.N=(int)width; // this is the wrong way round as I have row major ordering and imagemagick column major
    I.M=(int)height;
    I.I=malloc(I.N*I.M*sizeof(double));

    // extract the data
    MagickExportImagePixels(wand, 0, 0, I.N, I.M, "I", DoublePixel, I.I);
    
    TransposeFloatImage(&I);
    
    
    // wrap it up
	if(wand)
		wand = DestroyMagickWand(wand);	
	MagickWandTerminus();
	
	return I;
}

void FloatimageWrite(char *fn, image I, int norm, double min, double max)
{
	
	MagickWand *wand = NULL;
	PixelWand *p_wand = NULL;
	MagickBooleanType status;
	double *data;
	image *P;
	image N;
	
	// init lib
	MagickWandGenesis();
	wand=NewMagickWand();
	
	// init the new image
	p_wand = NewPixelWand();
	PixelSetColor(p_wand,"white");
	status=MagickNewImage(wand,I.M,I.N,  p_wand);
	
	if (status == MagickFalse)
		ThrowWandException(wand);
	
	if (norm)
	{
		N=NormalizeImageRange(I,min,max);
		P=&N;
	}
	else
		P=&I;
	
	// import the data
    data=TransposeImageData(*P);
	status=MagickImportImagePixels(wand,0,0,I.M,I.N,"I",DoublePixel,data);
	
	if (status == MagickFalse)
		ThrowWandException(wand);
		
	// write the image
	status=MagickWriteImage(wand,fn);
	
	// Diagnose any error
	if (status == MagickFalse)
		ThrowWandException(wand);
		
	// wrap it up
	if(wand)
		wand = DestroyMagickWand(wand);
	if (p_wand)
		p_wand = DestroyPixelWand(p_wand);

	free(data);
	if (norm)
		FreeImage(&N);
	MagickWandTerminus();
}

void FloatimageDisplay(image I, int norm, double min, double max)
{
	
	MagickWand *wand = NULL;
	PixelWand *p_wand = NULL;
	MagickBooleanType status;
	double *data;
	image *P;
	image N;
	
	// init lib
	MagickWandGenesis();
	wand=NewMagickWand();
	
	// init the new image
	p_wand = NewPixelWand();
	PixelSetColor(p_wand,"white");
	status=MagickNewImage(wand,I.M,I.N,  p_wand);
	
	if (status == MagickFalse)
		ThrowWandException(wand);
	
	if (norm)
	{
		N=NormalizeImageRange(I,min,max);
		P=&N;
	}
	else
		P=&I;
	
	// import the data
    data=TransposeImageData(*P);
	status=MagickImportImagePixels(wand,0,0,I.M,I.N,"I",DoublePixel,data);
	
	if (status == MagickFalse)
		ThrowWandException(wand);
		
	// display the image
	MagickDisplayImage(wand, ":0");
		
	// wrap it up
	if(wand)
		wand = DestroyMagickWand(wand);
	
    free(data);
	if (norm)
		FreeImage(&N);
	MagickWandTerminus();
}

/* dump image to stdout */
void Floatimage2stdout(image I)
{
	int i,j;
	for (i=0;i<I.N;i++)
	{
		for (j=0;j<I.M;j++)
			printf("%e\t",I.I[INDEX(i,j,I.N)]);
		printf("\n");
	}
	printf("\n");
}

/* dump image to plain text file */
void FloatimageTXTWrite(char *fn, image I)
{
	FILE *f;
	int i,j;
	
	if ((f=fopen(fn,"w"))==NULL)
	{
		
		// ERRORFLAG ERROFILEW  "cannot open file for writing"
		AddErr(ERROFILEW);
		return;
	}
	for (i=0;i<I.N;i++)
	{
		for (j=0;j<I.M;j++)
			fprintf(f,"%e\t",I.I[INDEX(i,j,I.N)]);
		fprintf(f,"\n");
	}
	fclose(f);
}

double *GetLine(FILE *f, int *M)
{
	double *Line;
	char *word;
	int Na=10;
	int i=0;
	char ch='a';
	
	Line=malloc(Na*sizeof(double));
	while ((feof(f)==0)&&(ch!='\n'))
	{
		/* get and allocate a word */
		if (fscanf(f,"%ms", &word)==1)
		{
			if (*word!='#')
			{
				/* put it in the data array */
				Line[i]=atof(word);
				i++;
				if (i==Na)
				{
					Na+=10;
					Line=realloc(Line,Na*sizeof(double));
				}	
			}
			/* free the word for the next run */
			free(word);
		}
		// clear whitespace but no endlines
		while ((feof(f)==0)&&(fscanf(f,"%[ \t]", &ch)==1));
		/* check for an endline */
		fscanf(f,"%[\n]",&ch);
	}
	Line=realloc(Line,i*sizeof(double));
	(*M)=i;
	return Line;
}


image FloatimageTXTRead(char *fn)
{
	FILE *f;
	image I = {NULL, 0, 0};
	double *Line;
	double **D;
	int M, Ml=0;
	int N=0, Na=10, i, j;
	
	if ((f=fopen(fn,"r"))==NULL)
	{
		// ERRORFLAG ERRIFILER  "cannot open txt file for reading"
		AddErr(ERRIFILER);
		return I;
	}
	do
	{
		Line=GetLine(f, &Ml);
		if (!Ml)
			free(Line);
	} while ((feof(f)==0)&&(Ml==0));
	
	if (Ml==0)
	{
		// ERRORFLAG ERRTXTINNODATA  "no usable data in file"
		AddErr(ERRTXTINNODATA);
		return I;
	}
	M=Ml;
	D=malloc(Na*sizeof(double *));
	D[N]=Line;
	N++;
	while ((feof(f)==0))
	{		
		Line=GetLine(f, &Ml);
		if (!Ml)
			free(Line);
		else if (Ml!=M)	
		{
			// ERRORFLAG ERRTXTINDIFFLL  "different line lengths in image file"
			AddErr(ERRTXTINDIFFLL);
			return I;
		}
		else
		{
			D[N]=Line;
			N++;
			if (N==Na)
			{
				Na+=10;
				D=realloc(D,Na*sizeof(double *));
			}			
		}
	}
	printf("loaded %dx%d image\n", N, M);
	I.M=M;
	I.N=N;	
    I.I=malloc(I.N*I.M*sizeof(double));	
	for (i=0;i<I.N;i++)
	{
		for (j=0;j<I.M;j++)
			I.I[INDEX(i,j,I.N)]=D[i][j];
		free(D[i]);
	}
	
	fclose(f);
	free(D);
	return I;
}

