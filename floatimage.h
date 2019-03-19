typedef struct image {
	double *I;
	int N, M;
} image;

image FloatimageRead(char *fn);
void FloatimageWrite(char *fn, image I, int norm, double min, double max);
void FloatimageDisplay(image I, int norm, double min, double max);
void FloatimageTXTWrite(char *fn, image I);
void PlotImage(image I);
void TransposeFloatImage(image I);
void FreeImage(image *I);
void EL2Vj(image I, double T);
void Vj2EL(image I, double T);
image NormalizeImageRange(image I, double min, double max);
