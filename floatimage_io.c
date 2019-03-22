#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <wand/MagickWand.h>
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
		fprintf(stderr, "Error: cannot open file %s for reading\n", fn);
		exit(1);
	}
	for (i=0;i<I.N;i++)
	{
		for (j=0;j<I.M;j++)
			fprintf(f,"%e\t",I.I[INDEX(i,j,I.N)]);
		fprintf(f,"\n");
	}
	fclose(f);
}
