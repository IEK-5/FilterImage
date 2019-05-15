#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "error.h"
#include "floatimage.h"

/* transpose image data: created a newly allocated image data transposed */
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

/* transpose an image */
void TransposeFloatImage(image *I)
{
	double *AA;
	int D;
	AA=TransposeImageData(*I);
	free(I->I);
	I->I=AA;
	D=I->M;
	I->M=I->N;
	I->N=D;
}

void FreeImage(image *I)
{
	I->N=0;
	I->M=0;
	free(I->I);
	I->I=NULL;
}

/* normalize image range to fit between min and max */
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
		// ERRORFLAG ERRNORMIM  "cannot normalize float image"
		AddErr(ERRNORMIM);
		f=10;
		c=0;
	}
	else
	{
		f=(max-min)/(dmax-dmin);
		c=min-f*dmin;
	}
	R.I=malloc(I.N*I.M*sizeof(double));
	R.N=I.N;
	R.M=I.M;
	for (i=0;i<I.N*I.M;i++)
		R.I[i]=f*I.I[i]+c;
	printf("scale: %e\t shift(0 @): %e)\n", f, c);
	return R;
}

/* convert EL intensity to junction voltage */
#define _kb 8.617330350e-5
void EL2Vj(image I, double T)
{
	int i;
	for (i=0;i<I.N*I.M;i++)
		I.I[i]=_kb*T*log(I.I[i]);
}
/* convert junction voltage to EL intensity */
void Vj2EL(image I, double T)
{
	int i;
	for (i=0;i<I.N*I.M;i++)
		I.I[i]=exp(I.I[i]/(_kb*T));		
}

/* returns a test image with I(x,y)=x^p+y^p */
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

image DupImage(image A)
{
	image R;
	int i;
	R.N=A.N;
	R.M=A.M;
	R.I=malloc(R.N*R.M*sizeof(double));
	for (i=0;i<R.N*R.M;i++)
		C.I[i]=A.I[i];	
}

void AddImages(image A, image B, double fa, double fb, image R)
{
	int i;
	if ((A.N!=B.N)||(A.M!=B.M)||(A.N!=C.N)||(A.M!=C.M))
	{
		// ERRORFLAG ERRADDIM  "cannot add images, dimensions do not match"
		AddErr(ERRADDIM);
		return;
	}
	for (i=0;i<C.N*C.M;i++)
		C.I[i]=fa*A.I[i]+fb*B.I[i];
}

void MultImages(image A, image B, double fa, double fb, image R)
{
	int i;
	if ((A.N!=B.N)||(A.M!=B.M)||(A.N!=C.N)||(A.M!=C.M))
	{
		// ERRORFLAG ERRMULTIM  "cannot multiply images, dimensions do not match"
		AddErr(ERRMULTIM);
		return;
	}
	for (i=0;i<C.N*C.M;i++)
		C.I[i]=fa*A.I[i]*fb*B.I[i];
}

void DivImages(image A, image B, double fa, double fb, image R)
{
	int i;
	if ((A.N!=B.N)||(A.M!=B.M)||(A.N!=C.N)||(A.M!=C.M))
	{
		// ERRORFLAG ERRDIVIM  "cannot divide images, dimensions do not match"
		AddErr(ERRDIVIM);
		return;
	}
	for (i=0;i<C.N*C.M;i++)
		C.I[i]=fa*A.I[i]/(fb*B.I[i]);
}
	
void ScaMultImage(image A, double fa)
{	
	int i;
	for (i=0;i<A.N*A.M;i++)
		A.I[i]*=fa;
}
void ScaDivImage(image A, double fa)
{	
	int i;
	for (i=0;i<A.N*A.M;i++)
		A.I[i]=fa/A.I[i];
}

void ScaAddImage(image A, double fa)
{	
	int i;
	for (i=0;i<A.N*A.M;i++)
		A.I[i]+=fa;
}
