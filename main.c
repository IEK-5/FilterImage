#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "floatimage.h"
#include "filter.h"

//#define DEBUG	
int main(int argc, char **argv)
{
	image Iin, Iout;
#ifdef DEBUG
	int N=3, m=3, m_nabla=2;// apply a second order nabla filter based on a 4-th order polynomal 
	Iin=TestImage(100,200,2);
	//Iout=FFT_PolynomalFilter(Iin, N, m, m_nabla, 'n'); // apply the filter
	Iout=PolynomalExtremaLocator(Iin, 10, 2);
	Floatimage2stdout(Iout);
#else	
	int N=10, m=4, m_nabla=2;// apply a second order nabla filter based on a 4-th order polynomal 
	
	if (argc!=3)
	{
		fprintf(stderr, "Usage: %s input-image output-image\n", argv[0]);
		exit(1);
	}
	Iin=FloatimageRead(argv[1]);
	
	EL2Vj(Iin, 300); // convert to junction voltages	
	//Floatimage2stdout(Iin);
	Iout=FFT_PolynomalExtremaLocator(Iin, 10, 7);
	Floatimage2stdout(Iout);
	//Iout=FFT_PolynomalFilter(Iin, N, m, m_nabla, 'n'); // apply the filter
	//FloatimageDisplay(Iin, 1,0,1);
	//FloatimageDisplay(Iout, 1,0,1);
	//FloatimageWrite(argv[2], Iout, 1,0,1);
#endif
	FreeImage(&Iin);
	FreeImage(&Iout);
	return 0;
}
