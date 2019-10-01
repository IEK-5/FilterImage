#ifndef FLOATIMAGE_IO_H
#define FLOATIMAGE_IO_H
image FloatimageRead(char *fn);
void FloatimageWrite(char *fn, image I, int norm, double min, double max);
void FloatimageDisplay(image I, int norm, double min, double max);
void FloatimageTXTWrite(char *fn, image I);
image FloatimageTXTRead(char *fn);
void Floatimage2stdout(image I);
#endif
