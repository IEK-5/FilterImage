/* FilterImage - A 2D Savitzky-Golay Image Filtering tool 
 * 
 * Copyright (C) 2021  Forschungszentrum Juelich GmbH
 * B. E. Pieters, E. Sovetkin, and  M. Gordon
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
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

extern double FILTER_EPS;

void PrintMat(double *A, int N, int M);
void FreeFilter(filter *F);
void FreeFilterSet(filterset *F);

// computing & handling filters
filter PartDeriv2DSum(int nn, int ns, int nw, int ne, int m, int *dmx, int *dmy, double *f, int N, char method);
void SaveFilter(char *fn, filter F);
filter LoadFilter(char *fn);

filterset DerivOperatorSet2D(int nn, int ns, int nw, int ne, int m, int *dmx, int *dmy, double *f, int N, char method);
void SaveFilterSet(char *fn, filterset F);
filterset LoadFilterSet(char *fn);

// filtering images using plain (slow!) convolution sums
image PL_ApplyFilter(image I, int stepy, int stepx,  filterset F);

// fft based filters, fast but with edge effects
image FFT_ApplyFilter(image I, int stepy, int stepx,  filter F);

// mix of plain and FFT routines, try to be fast without edge effects
image ApplyFilter(image I, int stepy, int stepx,  filterset F);

/************************** Compact Interface, without seperate filters *************************************/
image PL_PolynomalFilter(image I, int ny, int nx, int stepy, int stepx, int m, int *dmx, int *dmy, double *f, int Nd, char method);
image FFT_PolynomalFilter(image I, int ny, int nx, int stepy, int stepx, int m, int *dmx, int *dmy, double *f, int Nd, char method);
image PolynomalFilter(image I, int ny, int nx, int stepy, int stepx, int m, int *dmx, int *dmy, double *f, int Nd, char method);
#endif
