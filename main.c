#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "floatimage.h"
#include "filter.h"

	
int main(int argc, char **argv)
{
	image Iin, Iout;
	int N=10, m=4, m_nabla=2;
	if (argc!=3)
	{
		fprintf(stderr, "Usage: %s input-image output-image\n", argv[0]);
		exit(1);
	}
	Iin=FloatimageRead(argv[1]);
	
	EL2Vj(Iin, 300); // convert to junction voltages	
	Iout=PolynomalFilter(Iin, N, m, m_nabla, 'n');
	NormalizeImageRange(Iout, 0, 1);
	FloatimageWrite(argv[2], Iout);
	FreeImage(&Iin);
	FreeImage(&Iout);
	return 0;
}
