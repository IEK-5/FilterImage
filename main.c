#include <stdlib.h>
#include <stdio.h>
#include "parser.h"

int main(int argc, char **argv)
{
	ParseCommandline(argc-1, argv+1);
}
