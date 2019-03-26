#ifndef FLOATIMAGE_H
#define FLOATIMAGE_H
typedef struct image {
	double *I; // images must in column major format
	int N, M;
} image;

#define INDEX(R,C,N)  ((C)*(N)+(R)) // column major format
#define ROW(I,N)  ((I)%(N))
#define COL(J,N)  ((J)/(N))

double * TransposeImageData(image I);
void TransposeFloatImage(image *I);
void FreeImage(image *I);
void EL2Vj(image I, double T);
void Vj2EL(image I, double T);
image NormalizeImageRange(image I, double min, double max);
image TestImage(int N, int M, double p);
#endif