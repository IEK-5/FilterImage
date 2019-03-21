CFLAGS = -Wall -pedantic -Ofast -flto `MagickWand-config --cflags`
LFLAGS =  -lm -flto -llapack -lopenblas -lfftw3 `MagickWand-config --ldflags --libs`
CC = gcc
SRC=floatimage.c filter.c main.c
OBJ=floatimage.o filter.o main.o
HDR=floatimage.h filter.h
TARGET=FilterImage

$(TARGET): $(OBJ)
	$(CC) -o $(TARGET) $(OBJ) $(LFLAGS) 
filter.o: floatimage.h filter.h
main.o: floatimage.h filter.h
clean:
	rm -f *.o $(TARGET)
