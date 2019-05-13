CFLAGS = -Wall -pedantic -Ofast -flto -fPIC `MagickWand-config --cflags`
LFLAGS =  -lm -flto -llapack -lopenblas -lfftw3 `MagickWand-config --ldflags --libs`
CC = gcc
SRC=floatimage.c floatimage_io.c filter.c error.c main.c parser.c variables.c
OBJ=floatimage.o floatimage_io.o filter.o main.o error.o parser.o variables.o
HDR=floatimage.h floatimage_io.h filter.h error.h parser.h variables.h
TARGET=FilterImage

$(TARGET): $(OBJ)
	$(CC) -o $(TARGET) $(OBJ) $(LFLAGS) 
filterimage.oct: filterimage.cc floatimage.o filter.o error.o
	mkoctfile -v  filterimage.cc  floatimage.o filter.o error.o
filter.o: floatimage.h filter.h error.h errorflags.h errormessages.h
floatimage_io.o: floatimage.h filter.h error.h errorflags.h errormessages.h
floatimage.o: error.h errorflags.h errormessages.h
main.o: parser.h
error.o: errorflags.h errormessages.h
variables.o: floatimage.h filter.h error.h
parser.o: parsedefs.h floatimage_io.h floatimage.h filter.h error.h variables.h
errorflags.h errormessages.h: gen_errorflags.sh filter.c floatimage.c floatimage_io.c error.c parser.c variables.c
	${SHELL} gen_errorflags.sh filter.c floatimage.c floatimage_io.c error.c parser.c variables.c
parsedefs.h: gen_parserflags.sh parser.c
	${SHELL} gen_parserflags.sh parser.c
clean:
	rm -f *.o $(TARGET) errorflags.h errormessages.h parserdef.h
