#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <wand/MagickWand.h>
#include "floatimage.h"

double * TransposeImageData(image I)
{
	double *AA;
	int i, j;
	AA=malloc(I.N*I.M*sizeof(double));
	
	for (i=0;i<I.N;i++)
		for (j=0;j<I.M;j++)
			AA[INDEX(j,i,I.M)]=I.I[INDEX(i,j,I.N)];
	return AA;
}

void TransposeFloatImage(image *I)
{
	double *AA;
	int i, D;
	AA=TransposeImageData(*I);
	for (i=0;i<I->N*I->M;i++)
		I->I[i]=AA[i];
	D=I->M;
	I->M=I->N;
	I->N=D;
	free(AA);
}

void FreeImage(image *I)
{
	I->N=0;
	I->M=0;
	free(I->I);
	I->I=NULL;
}

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

image NormalizeImageRange(image I, double min, double max)
{
	double dmin, dmax, f, c;
	image R;
	int i;
	dmin=dmax=I.I[0];
	for (i=1;i<I.N*I.M;i++)
	{
		if (I.I[i]<dmin)
			dmin=I.I[i];
		if (I.I[i]>dmax)
			dmax=I.I[i];
	}
	if (fabs(dmax-dmin)/fabs(dmin+dmax)<1e-10)
	{
		fprintf(stderr, "Error: cannot normalize float image %e %e\n", dmin, dmax);
		exit(1);
	}
	f=(max-min)/(dmax-dmin);
	c=min-f*dmin;
	R.I=malloc(I.N*I.M*sizeof(double));
	R.N=I.N;
	R.M=I.M;
	for (i=0;i<I.N*I.M;i++)
		R.I[i]=f*I.I[i]+c;
	printf("scale: %e\t shift(0 @): %e)\n", f, c);
	return R;
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
#define _kb 8.617330350e-5
void EL2Vj(image I, double T)
{
	int i;
	for (i=0;i<I.N*I.M;i++)
		I.I[i]=_kb*T*log(I.I[i]);
}
void Vj2EL(image I, double T)
{
	int i;
	for (i=0;i<I.N*I.M;i++)
		I.I[i]=exp(I.I[i]/(_kb*T));		
}

image TestImage(int N, int M, double p)
{
	int i, j;
	image R;
	R.N=N;
	R.M=M;
	R.I=malloc(N*M*sizeof(double));	
	for (i=0;i<N;i++)
		for (j=0;j<M;j++)
			R.I[INDEX(i,j,N)]=pow((double)(i-N/2),p)+pow((double)(j-M/2),p);
	return R;	
} 

/*
int main(int argc, char **argv)
{
	image I;
	if (argc!=3)
	{
		fprintf(stderr, "xx in out\n");
		exit(1);
	}
	I=MagickRead(argv[1]);
	PlotImage(I);
	MagickWrite(argv[2],I);
	return(0);
}*/

