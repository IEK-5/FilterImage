CFLAGS = -Wall -pedantic -Ofast -flto `MagickWand-config --cflags`
LFLAGS =  -lm -flto -llapack -lopenblas -lfftw3 `MagickWand-config --ldflags --libs`
CC = gcc
SRC=floatimage.c floatimage_io.c filter.c main.c
OBJ=floatimage.o floatimage_io.o filter.o main.o
HDR=floatimage.h floatimage_io.h filter.h
TARGET=FilterImage

$(TARGET): $(OBJ)
	$(CC) -o $(TARGET) $(OBJ) $(LFLAGS)
filter.o: floatimage.h filter.h
floatimage_io.o: floatimage.h filter.h
main.o: floatimage_io.h floatimage.h filter.h
clean:
	rm -f *.o $(TARGET)
