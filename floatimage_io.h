image FloatimageRead(char *fn);
void FloatimageWrite(char *fn, image I, int norm, double min, double max);
void FloatimageDisplay(image I, int norm, double min, double max);
void FloatimageTXTWrite(char *fn, image I);
void Floatimage2stdout(image I);