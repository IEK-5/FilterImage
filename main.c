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
