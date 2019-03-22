#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "floatimage.h"
#include "filter.h"

void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *A, int *LDA, double *B, int *LDB, double *BETA, double *C,int *LDC); 

void Transpose(double *A, int N, int M)
{
	double *AA;
	int i, j;
	AA=malloc(N*M*sizeof(double));
	
	for (i=0;i<N;i++)
		for (j=0;j<M;j++)
			AA[INDEX(j,i,M)]=A[INDEX(i,j,N)];
	for (i=0;i<N*M;i++)
		A[i]=AA[i];
	free(AA);
}

void MMultABT(int Na, int Ma, int Nb, int Mb, double *A, double *B, double *R)
{
	char transA = 'N', transB = 'T';
	double one = 1.0, zero = 0.0;
	if (Ma!=Mb)
	{
		fprintf(stderr,"Error: cannot multiply matrices, domensions do not match\n");
		exit(1);
	}

	dgemm_(&transA, &transB, &Na, &Nb, &Ma, &one, A,&Na, B, &Nb, &zero, R, &Na);
}

void ATA_(int N, int M, double *a, double *r)
/* at*a=r */
{
	int i,j, k, z;	
	for (i=0;i<M;i++)
	{
		for (j=0;j<M;j++)
		{
			z=INDEX(i,j,M);
			r[z]=0;
			for (k=0;k<N;k++)
				r[z]+=a[INDEX(k,i,N)]*a[INDEX(k,j,N)];
		}
	}
}

void ATA(int N, int M, double *a, double *r)
/* at*a=r */
{
	char transA = 'T', transB = 'N';
	double one = 1.0, zero = 0.0;

	dgemm_(&transA, &transB, &M, &M, &N, &one, a,&N, a, &N, &zero, r, &M);
}

void PrintMat(double *A, int N, int M)
{
	int i,j;
	printf("\n");
	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
			printf("%e\t",A[INDEX(i,j,N)]);
		printf("\n");
	}
	printf("\n");
}

void inverse(double* A, int N)
{
    int *IPIV;
    int LWORK = N*N;
    double *WORK;
    int INFO;
    
    IPIV = malloc((N+1)*sizeof(int));
    WORK = malloc(LWORK*sizeof(double));

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    free(IPIV);
    free(WORK);
}

double *LLS(int nn, int ns, int nw, int ne, int m)
/* 
 * nn: number values to the north
 * ns: number values to the south
 * nw: number values to the west
 * ne: number values to the east
 * m: order of the polynomal
 * Returns the linear least squares soltution matrix : (A^TA)^-1 A^T
 */

{
	double *A, *AA, *a;
	int N, M;
	int I,J,i,j;
	if (nn<0)
	{
		fprintf(stderr,"Error: number of north elements may not be negative\n");
		exit(1);
	}
	if (ns<0)
	{
		fprintf(stderr,"Error: number of south elements may not be negative\n");
		exit(1);
	}
	if (nw<0)
	{
		fprintf(stderr,"Error: number of west elements may not be negative\n");
		exit(1);
	}
	if (ne<0)
	{
		fprintf(stderr,"Error: number of east elements may not be negative\n");
		exit(1);
	}
	
	N=(nn+ns+1)*(nw+ne+1);	
	M=((nn+ns>0)+(ne+nw>0))*m+1;
	
	if (M>N)
	{
		fprintf(stderr,"Error: number of datapoints not sufficient for the polynomal order\n");
		exit(1);
	}
	
	A=malloc((N*M)*sizeof(double));
	AA=malloc((M*M)*sizeof(double));
	i=0;
	for (I=-nn;I<=ns;I++)
	{
		for (J=-nw;J<=ne;J++)
		{
			int k=0;
			A[INDEX(i,k++,N)]=1.0;
			if (nn+ns>0)
				for (j=1;j<m+1;j++)
					A[INDEX(i,k++,N)]=pow((double)I,(double)j);
			if (ne+nw>0)
				for (j=1;j<m+1;j++)
					A[INDEX(i,k++,N)]=pow((double)J,(double)j);
			i++;
		}
	}
	ATA(N, M, A, AA); // AA=A^T A
	inverse(AA, M); // AA=(A^T A)^-1
	a=malloc((N*M)*sizeof(double));
	
	MMultABT(M, M, N, M, AA, A, a); // a=(A^T A)^-1 A^T
	
	free(A);
	free(AA);
	return a;
}


filter PartDeriv2D(int nn, int ns, int nw, int ne, int deriv_m, int m, double fx, double fy)
/* 
 * nn: number values to the north
 * ns: number values to the south
 * nw: number values to the west
 * ne: number values to the east
 * deriv_m: Order of the differentaitor operator (may also be 0, which results in smoothing)
 * m: order of the polynomal
 * d: direction, can be 'x' (partial derivative to x), 'y' (partial derivative to y), 'n' (nabla operator, sum of partial derivatives to x and y)
 * returns the filter coefficients in column major format (vector of size (nn+ns*1)*(nw+ne*1))
 */
{
	
	filter F;
	double *a;
	double dd;
	int N, M;
	int i;

	a=LLS(nn, ns, nw, ne, m);
	
	N=(nn+ns+1)*(nw+ne+1);	
	M=((nn+ns>0)+(ne+nw>0))*m+1;
	
	
	dd=1.0;
	for (i=0;i<deriv_m;i++)
		dd*=(double)(i+1);
	
	F.nn=nn;
	F.ns=ns;
	F.ne=ne;
	F.nw=nw;	
	F.F=calloc(N,sizeof(double));
	
#define EPS 1e-10	
	if (fabs(fx)+fabs(fy)<EPS)
	{
		fprintf(stderr,"Error: you have to specify a non zero length direction vector\n");
		exit(1);
	}
	
	if (fabs(fx)/(fabs(fx)+fabs(fy)) > EPS)
	{
		if ((nw+ne)<=0)
		{
			fprintf(stderr,"Error: cannot make derivative in X direction without east or west elements\n");
			exit(1);
		}
		if (nn+ns>0)
		{
			for (i=0;i<N;i++)
				F.F[i]+=fx*dd*a[INDEX(m+deriv_m,i,M)];
		}
		else
		{
			for (i=0;i<N;i++)
				F.F[i]+=fx*dd*a[INDEX(deriv_m,i,M)];
		}
	}
	if (fabs(fy)/(fabs(fx)+fabs(fy)) > EPS)
	{
		if ((nn+ns)<=0)
		{
			fprintf(stderr,"Error: cannot make derivative in Y direction without north or south elements\n");
			exit(1);
		}
		for (i=0;i<N;i++)
			F.F[i]+=fy*dd*a[INDEX(deriv_m,i,M)];
	}	
	Transpose(F.F, (ne+nw+1), (nn+ns+1)); // I suppose I could modify LLS to avoid this step...
	free(a);
	return F;
}

void PrintFilter(filter F)
// print filter coefficients
{
	printf("x:-%d ... %d y: -%d ... %d\n", F.nw, F.ne, F.nn, F.ns);
	PrintMat(F.F, (F.nn+F.ns+1), (F.nw+F.ne+1));
}

void FreeFilter(filter *F)
{
	F->nn=0;
	F->ns=0;
	F->ne=0;
	F->nw=0;
	free(F->F);
	F->F=NULL;
}

filterset DerivOperatorSet2D(int nn, int ns, int nw, int ne, int deriv_m, int m, double fx, double fy)
/* 
 * nn: number values to the north
 * ns: number values to the south
 * nw: number values to the west
 * ne: number values to the east
 * deriv_m: Order of the differentaitor operator (may also be 0, which results in smoothing)
 * m: order of the polynomal
 * d: direction, can be 'x' (partial derivative to x), 'y' (partial derivative to y), 'n' (nabla operator, sum of partial derivatives to x and y)
 * Computes a set of filters. Any single filter will lead to problems at the edges of an image
 * The computed filter set allows filtering the image up to all edges of the image
 */
{
	int i,j;
	int n, s, w, e;
	int N, M;
	filterset F;
	N=(nn+ns+1);
	M=(nw+ne+1);
	
	F.nn=nn;
	F.ns=ns;
	F.ne=ne;
	F.nw=nw;
	F.set=malloc(N*M*sizeof(filter));
	for (i=-nn;i<=ns;i++)
		for (j=-nw;j<=ne;j++)
		{		
			n=nn;
			s=ns;
			e=ne;
			w=nw;
			if (i<0)
				n+=i;
			if (i>0)
				s-=i;
			if (j<0)
				w+=j;
			if (j>0)
				e-=j;
			F.set[(nn+i)*M+(nw+j)]=PartDeriv2D(n, s, w, e, deriv_m, m, fx, fy);
			//if ((n==nn)&&(s==ns)&&(e==ne)&&(w==nw))
			//	PrintMat(F.set[(nn+i)*M+(nw+j)].F, (n+s+1), w+e+1);
		}
	return F;	
}


void FreeFilterSet(filterset *F)
{
	int i,j;
	for (i=-F->nn;i<=F->ns;i++)
		for (j=-F->nw;j<=F->ne;j++)
			FreeFilter(F->set+(F->nn+i)*(F->nw+F->ne+1)+(F->nw+j));
	F->nn=0;
	F->ns=0;
	F->ne=0;
	F->nw=0;
	free(F->set);
	F->set=NULL;
}


void PL_ApplyFilterRange(image I, image R, int row0, int rown, int col0, int colm, int stepy, int stepx, filterset F)
{
	int i, j, k,l;
	for (i=row0;i<=rown;i++)
		for (j=col0;j<=colm;j++)
		{
			double s=0;
			filter f;
			int ii=0, jj=0;
			/* OK this part is a bit of a disaster
			 * We need to select the right filter such that we do not get out of the image bounds
			 */
			if (i-F.nn*stepy<0)
				ii=i-F.nn*stepy-i%stepy;
			if (i+F.ns*stepy-I.N+1>0)
				ii=i+F.ns*stepy-I.N+1+((I.N-1-i)%stepy);
				
			if (j-F.nw*stepx<0)
				jj=j-F.nw*stepx-(j%stepx);
			if (j+F.ne*stepx-I.M+1>0)
				jj=j+F.ne*stepx-I.M+1+((I.M-1-j)%stepx);
				
			f=F.set[(F.nn+ii/stepy)*(F.nw+F.ne+1)+(F.nw+jj/stepx)];
			
			for (k=-f.nn;k<=f.ns;k++)
				for (l=-f.nw;l<=f.ne;l++)
					s+=f.F[INDEX(k+f.nn,l+f.nw,f.nn+f.ns+1)]*I.I[INDEX(i+k*stepy,j+l*stepx,I.N)];
			R.I[INDEX(i,j,I.N)]=s;
		}
}




image PL_ApplyFilter(image I, int stepy, int stepx, filterset F)
/* apply a filter set using a plain convolution */
{
	image R;
	R.I=malloc(I.N*I.M*sizeof(double));
	R.N=I.N;
	R.M=I.M;	
	PL_ApplyFilterRange(I, R, 0, I.N-1, 0, I.M-1, stepy, stepx, F);
	return R;
}

image PL_PolynomalFilter(image I, int ny, int nx, int stepy, int stepx, int m, int deriv_m, double fx, double fy)
/* the plain convolution filter
 * input:
 * I: 		image in column major format 
 * N,M: 	size of the image
 * nsurr:	number of elements to use in north, south, west, and east direction
 * deriv_m:	which derivative to compute
 * d:		dpecification of direction, x, y, or n (x+y)
 * Returns a filtered image.
 */
{	
	filterset F;
	image R;
	F=DerivOperatorSet2D(ny, ny, nx, ny, deriv_m, m, fx, fy);
	R=PL_ApplyFilter(I, stepy, stepx, F);
	FreeFilterSet(&F);
	return R;
}




image FFT_ApplyFilter(image I, int stepy, int stepx, filter F)
/* apply a filter set using the FFT for a fast convolution */
{
    fftw_complex *FI=NULL, *FP=NULL, *P;
	fftw_plan plan_inverse;
	fftw_plan plan_forward;
    double scale=1.0/((double)(I.M*I.N));
    double *PSF, a[2];
    image R;
	int i, j, ii,jj, ww;
	// create a filter image
	PSF=calloc(I.N*I.M, sizeof(double));
	for(i=-F.nn*stepy;i<=F.ns*stepx;i+=stepy)
		for (j=-F.nw*stepx;j<=F.ne*stepx;j+=stepx)
		{
			ii=i;
			jj=j;
			if (i<0)
				ii+=I.N;
			if (j<0)
				jj+=I.M;			
			PSF[INDEX(ii,jj,I.N)]=F.F[INDEX(F.nn+i/stepy, F.nw+j/stepx,F.nn+F.ns+1)];
		}
	//PrintMat(PSF, I.N, I.M);
	//PrintMat(F.F, F.nn*F.ns+1, F.nw*F.ne+1);
    
    // allocate working space
	ww=2*(I.N/2+1);
	FP = fftw_malloc ( sizeof ( double ) * I.M * ww );
	FI = fftw_malloc ( sizeof ( double ) * I.M * ww );
    plan_forward=fftw_plan_dft_r2c_2d ( I.M, I.N, PSF, FP, FFTW_ESTIMATE );
    //plan_forward_i=fftw_plan_dft_r2c_2d ( I.M, I.N, I.I,   FI, FFTW_ESTIMATE );
    
    // transform
	fftw_execute(plan_forward);
	fftw_execute_dft_r2c(plan_forward, I.I, FI);
	// convolute
	P=malloc(I.M*ww*sizeof(double));
	for (j=0;j<I.M*(I.N/2+1); j++)
	{		
		a[0]=(FI[j][0]*FP[j][0]-FI[j][1]*FP[j][1]);
		a[1]=(FI[j][0]*FP[j][1]+FI[j][1]*FP[j][0]);
		P[j][0]=scale*a[0];
		P[j][1]=scale*a[1];
	}
	// back transform
	R.N=I.N;
	R.M=I.M;
	R.I=malloc(I.N*I.M*sizeof(double));
	plan_inverse=fftw_plan_dft_c2r_2d ( I.M, I.N, P, R.I, FFTW_ESTIMATE );
	fftw_execute(plan_inverse);
	
	//cleanup
    free(FI);
    free(FP);
    free(P);
    free(PSF);
	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_inverse);
    
    return R;
}

image FFT_PolynomalFilter(image I, int ny, int nx, int stepy, int stepx, int m, int deriv_m, double fx, double fy)
/* same as above only now using an fft. 
 * downside: edge effects as the FFT treats the image as periodic
 * upside: This routine's computation time is independant on the filter size and thus this routine is faster for larger filters
 *         (for small filters, ~3x3, computation times are very similar)
 */
{	
    image R;
	filter F;
	F=PartDeriv2D(ny, ny, nx, nx, deriv_m, m, fx, fy);
	R=FFT_ApplyFilter(I, stepy, stepx, F);
	FreeFilter(&F);
    return R;
}


void PL_ApplyFilterEdge(image I, image R, int stepy, int stepx, filterset F)
/* apply a plain convolution filter to the edge */
{
	filter f0;
	f0=F.set[F.nn*(F.nw+F.ne+1)+F.nw];
	int north, south, west, east, Rlast, Clast;
	
	Rlast=I.N-1;
	Clast=I.M-1;
	
	north=f0.nn*stepy;
	south=Rlast-f0.ns*stepy;
	west=f0.nw*stepx;
	east=Clast-f0.ne*stepy,
	/* north edge */
	PL_ApplyFilterRange(I, R, 0      , north  , 0   , Clast, stepy, stepx, F);
	/* south edge */
	PL_ApplyFilterRange(I, R, south  , Rlast  , 0   , Clast, stepy, stepx, F);
	/* west edge */
	PL_ApplyFilterRange(I, R, north+1, south-1, 0   , west , stepy, stepx, F);
	/* east edge */
	PL_ApplyFilterRange(I, R, north+1, south-1, east, Clast, stepy, stepx, F);	
}


image ApplyFilter(image I, int stepy, int stepx, filterset F)
/* apply a filter set using a plain convolution */
{
	image R;
	R.I=malloc(I.N*I.M*sizeof(double));
	R.N=I.N;
	R.M=I.M;	
	/* first do an FFT convolution */
	R=FFT_ApplyFilter(I, stepy, stepx,F.set[F.nn*(F.nw+F.ne+1)+F.nw]);	
	/* redo the edge with more care */
	PL_ApplyFilterEdge(I, R, stepy, stepx, F);
	return R;
}

image PolynomalFilter(image I, int ny, int nx, int stepy, int stepx, int m, int deriv_m, double fx, double fy)
/* use FFT for the center part of the image and the plain method for the edge
 * input:
 * I: 		image in column major format 
 * N,M: 	size of the image
 * nsurr:	number of elements to use in north, south, west, and east direction
 * deriv_m:	which derivative to compute
 * d:		dpecification of direction, x, y, or n (x+y)
 * Returns a filtered image.
 */
{	
	filterset F;
	image R;
	F=DerivOperatorSet2D(ny, ny, nx, ny, deriv_m, m, fx, fy);
	R=ApplyFilter(I, stepy, stepx, F);
	FreeFilterSet(&F);
	return R;
}
