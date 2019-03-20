
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
filter PartDeriv2D(int nn, int ns, int nw, int ne, int deriv_m, int m, char d);
filterset DerivOperatorSet2D(int nn, int ns, int nw, int ne, int deriv_m, int m, char d);

void PrintFilter(filter F);
void FreeFilter(filter *F);
void FreeFilterSet(filterset *F);


// filtering images
image ApplyFilter(image I, filterset F);
image PolynomalFilter(image I, int nsurr, int m, int deriv_m, char d);
image PolynomalExtremaLocator(image I, int nsurr, double r);


// fft based filters, faster for large filters but with edge effects
image FFT_ApplyFilter(image I, filter F);
image FFT_PolynomalFilter(image I, int nsurr, int m, int deriv_m, char d);
image FFT_PolynomalExtremaLocator(image I, int nsurr, double r);


