
typedef struct filter {
	int nn, ns, ne, nw;
	double * F;
} filter;

typedef struct filterset {
	int nn, ns, ne, nw;
	filter *set;
} filterset;
void PrintMat(double *A, int N, int M);

filter PartDeriv2D(int nn, int ns, int nw, int ne, int deriv_m, int m, char d);
filterset DerivOperatorSet2D(int nn, int ns, int nw, int ne, int deriv_m, int m, char d);

image ApplyFilter(image I, filterset F);
image PolynomalFilter(image I, int nsurr, int m, int deriv_m, char d);
image FFT_PolynomalFilter(image I, int nsurr, int m, int deriv_m, char d);

void PrintFilter(filter F);

void FreeFilter(filter *F);
void FreeFilterSet(filterset *F);


