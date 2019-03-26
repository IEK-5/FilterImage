#ifndef FILTER_H
#define FILTER_H
#include "floatimage.h"
typedef struct filter {
	int nn, ns, ne, nw;
	double * F;
} filter;

typedef struct filterset {
	int nn, ns, ne, nw;
	filter *set;
} filterset;
void PrintMat(double *A, int N, int M);

// computing & handling filters
filter PartDeriv2D(int nn, int ns, int nw, int ne, int deriv_m, int m, double fx, double fy);
filterset DerivOperatorSet2D(int nn, int ns, int nw, int ne, int deriv_m, int m, double fx, double fy);

// filtering images using plain (slow!) convolution sums
image PL_ApplyFilter(image I, int stepy, int stepx,  filterset F);
image PL_PolynomalFilter(image I, int ny, int nx, int stepy, int stepx, int m, int deriv_m, double fx, double fy);

// fft based filters, fast but with edge effects
image FFT_ApplyFilter(image I, int stepy, int stepx,  filter F);
image FFT_PolynomalFilter(image I, int ny, int nx, int stepy, int stepx, int m, int deriv_m, double fx, double fy);

// mix of plain and FFT routines, try to be fast without edge effects
image ApplyFilter(image I, int stepy, int stepx,  filterset F);
image PolynomalFilter(image I, int ny, int nx, int stepy, int stepx, int m, int deriv_m, double fx, double fy);
#endif

