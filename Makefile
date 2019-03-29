CFLAGS = -Wall -pedantic -Ofast -flto `MagickWand-config --cflags`
LFLAGS =  -lm -flto -llapack -lopenblas -lfftw3 `MagickWand-config --ldflags --libs`
CC = gcc
SRC=floatimage.c floatimage_io.c filter.c error.c main.c
OBJ=floatimage.o floatimage_io.o filter.o main.o error.o
HDR=floatimage.h floatimage_io.h filter.h error.h 
TARGET=FilterImage

$(TARGET): $(OBJ)
	$(CC) -o $(TARGET) $(OBJ) $(LFLAGS) 
filter.o: floatimage.h filter.h error.h errorflags.h errormessages.h
floatimage_io.o: floatimage.h filter.h error.h errorflags.h errormessages.h
floatimage.o: error.h errorflags.h errormessages.h
main.o: floatimage_io.h floatimage.h filter.h error.h
error.o: errorflags.h errormessages.h
errorflags.h errormessages.h: gen_errorflags.sh filter.c floatimage.c floatimage_io.c error.c
	${SHELL} gen_errorflags.sh filter.c floatimage.c floatimage_io.c error.c
clean:
	rm -f *.o $(TARGET) errorflags.h errormessages.h
