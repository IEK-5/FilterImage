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
image DupImage(image A);
void AddImages(image A, image B, double fa, double fb, image R);
void MultImages(image A, image B, double fa, double fb, image R);
void DivImages(image A, image B, double fa, double fb, image R);	
void ScaMultImage(image A, double fa);
void ScaDivImage(image A, double fa);
void ScaAddImage(image A, double fa);
void AbsImage(image A);
void ImageClip(image A, double th, int p);
#endif
