#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "parser.h"
#include "readlineshell.h"
#include "filter.h"
#include "error.h"
#include "variables.h"

#define MAXLINELEN 4096

void StripEndline(char *line)
{
	int i;
	i=strlen(line)-1;
	while ((i>0)&&((line[i]=='\n')||(line[i]=='\r')))
	{
		line[i]='\0';
		i--;
	}
}
void ParseFile(char *fn)
{
	FILE *f;
	char *line;	
	if ((f=fopen(fn,"r"))==NULL)
	{
		
		// ERRORFLAG ERRPFILEREAD  "cannot open input file for reading"
		AddErr(ERRPFILEREAD);
		return;
	}
	line=malloc(MAXLINELEN*sizeof(char));
    fgets(line, MAXLINELEN-1, f);
	while(feof(f)==0)
	{
		StripEndline(line);
		if (ParseComm(line))
			break;
		fgets(line, MAXLINELEN-1, f);
	}
	free(line);
}
int main(int argc, char **argv)
{
	int i;
	InitVars();
	for (i=1;i<argc;i++)
	{
		if (argv[i][0]=='-')
		{
			switch (argv[i][1])
			{
				case 'f':
					i++;
					if (i<argc)
						ParseFile(argv[i]);
					else
						fprintf(stderr,"Error: missing file after -f option\n");
					break;
				case 'i':
					shell();
					break;
				default:
					fprintf(stderr,"Error: unknown option %s\n", argv[i]);
			}
		}
		else
		
		if (ParseComm(argv[i]))
			break;
	}
	ClearVars();

	exit(EXITSTATE);
}
