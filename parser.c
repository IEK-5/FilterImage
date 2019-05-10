#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
/* local includes */
#include "parsedef.h"
#include "floatimage.h"
#include "floatimage_io.h"
#include "filter.h"
#include "error.h"

char * GetWord(char *in, char *word)
/* take string in and allocated string word
 * copy the first word of strin in to word
 * return a pointer to the rest of the string
 */
{
	char *end=in;
	char *out=word;
	int go=1;
	while (*end && go)
	{
		if (isblank(*end))
		{
			go=0;
		}
		else 
		{
			*out=*end;
			out++;
		}
		end++;
	}
	*out='\0';
	return end;
}

ParserFun LookupComm(char *key)
{
	int i=0;
	int lk;
	lk=strlen(key);
	while (KeyTable[i].key)
	{
		if (strlen(KeyTable[i].key)==lk)
			if (strncmp(KeyTable[i].key, key, lk)==0)
				return (KeyTable[i].fun);
		i++;
	}
	return NULL;
}


void ParseComm(char *in)
{
	char *key;
	char *arg;
	ParserFun fun;
	key=malloc(strlen(in)*sizeof(char));
	arg=GetWord(in, key);
	fun=LookupComm(key);
	if (fun)
		fun(arg);
	else
	{		
		// ERRORFLAG ERRCOMMUND  "command not defined"
		AddErr(ERRCOMMUND);
	}
	if (ERRORSTATE)
	{
		E_Messages();
		printf("%d-->%s\n", ERRORSTATE, in);
		ResetErrors();
	}
	free(key);
}

void GetOption(char *in, char *opt, char *word)
{
	char *start;
	char *opti;
	int len;
	len=(strlen(opt)+2);
	opti=malloc(len*sizeof(char));
	snprintf(opti,len,"%s=",opt);
	start=strstr(in, opti);
	if (!start)
	{
		fprintf(stderr,"Value  %s not found\n", opt);
		// ERRORFLAG ERROPTIONFAILS  "Required value not found in input line"
		AddErr(ERROPTIONFAILS);
		free(opti);
		return;
	}	
	GetWord(start+len-1, word);
	free(opti);	
}


/* global variables */	
image Iin =  {NULL, 0, 0};
image Iout = {NULL, 0, 0};


/* parse routines
 * for each parse routine add a parseflag:
 * '// PARSEFLAG <keyword> <function-name>'
 * Note: you can create aliases by adding several lines like this with 
 * different keywords but the same function.
 * a parse routine must be of type void and take one character string as input
 * The gen_parseflags will take care of all the definitions in parsedef.h
 */
 
// PARSEFLAG read ReadFile
// PARSEFLAG load ReadFile
void ReadFile(char *in)
{
	int go=1;
	while (*in && go)
	{
		if (!isblank(*in))
		{
			go=0;
		}
		else
			in++;
	}
	if (Iin.I)
		FreeImage(&Iin);	
	Iin=FloatimageRead(in);
}

// PARSEFLAG save SaveFile
void SaveFile(char *in)
{
	int go=1;
	if (!Iout.I)
	{
		// ERRORFLAG ERRNOOIMAGE  "cannot save output image, no output image computed"
		AddErr(ERRNOOIMAGE);
		return;
	}
	while (*in && go)
	{
		if (!isblank(*in))
		{
			go=0;
		}
		else
			in++;
	}	
	FloatimageWrite(in, Iout, 1,0,1);
}

// PARSEFLAG txtsave TXTSaveFile
void TXTSaveFile(char *in)
{
	if (!Iout.I)
	{
		AddErr(ERRNOOIMAGE);
		return;
	}
	if (in)
	{
		if (strncmp(in, "-", 2)==0)
			Floatimage2stdout(Iout);
		else
			FloatimageTXTWrite(in, Iout);
	}
}

// PARSEFLAG ELtoVj ELtoVj
void ELtoVj(char *in)
{
	char *arg;
	char *temp;
	if (!Iin.I)
	{
		// ERRORFLAG ERRNOIIMAGE  "Cannot process input image, no image loaded"
		AddErr(ERRNOIIMAGE);
		return;
	}
	if (in)
	{
		temp=malloc(strlen(in)*sizeof(char));
		arg=GetWord(in, temp);
		EL2Vj(Iin, atof(temp)); // convert to junction voltages
		free(temp);
	}
}


// PARSEFLAG plain_filter PlainFilter
void PlainFilter(char *in)
{
	char *arg;
	char *opt;
	char *word;
	int Nx, Ny, stepy, stepx, m, dm;
	double fx, fy;
	
	if (!Iin.I)
	{
		AddErr(ERRNOIIMAGE);
		return;
	}
	
	// Parse options:
	word=malloc(strlen(in)*sizeof(char));
	
	GetOption(in, "Nx", word);
	Nx=atoi(word);	
	
	GetOption(in, "Ny", word);
	Ny=atoi(word);
	
	GetOption(in, "stepx", word);
	stepx=atoi(word);
	
	GetOption(in, "stepy", word);
	stepy=atoi(word);
	
	GetOption(in, "m", word);
	m=atoi(word);
	
	GetOption(in, "dm", word);
	dm=atoi(word);
	
	GetOption(in, "fx", word);
	fx=atof(word);
	
	GetOption(in, "fy", word);
	fy=atof(word);
	
	if (Iout.I)
		FreeImage(&Iout);
	Iout=PL_PolynomalFilter(Iin, Ny, Nx, stepy, stepx, m, dm, fx, fy);
}	
	
// PARSEFLAG filter Filter
void Filter(char *in)
{
	char *arg;
	char *opt;
	char *word;
	int Nx, Ny, stepy, stepx, m, dm;
	double fx, fy;
	
	if (!Iin.I)
	{
		AddErr(ERRNOIIMAGE);
		return;
	}
	
	// Parse options:
	word=malloc(strlen(in)*sizeof(char));
	
	GetOption(in, "Nx", word);
	Nx=atoi(word);	
	
	GetOption(in, "Ny", word);
	Ny=atoi(word);
	
	GetOption(in, "stepx", word);
	stepx=atoi(word);
	
	GetOption(in, "stepy", word);
	stepy=atoi(word);
	
	GetOption(in, "m", word);
	m=atoi(word);
	
	GetOption(in, "dm", word);
	dm=atoi(word);
	
	GetOption(in, "fx", word);
	fx=atof(word);
	
	GetOption(in, "fy", word);
	fy=atof(word);
	
	if (Iout.I)
		FreeImage(&Iout);
	Iout=PolynomalFilter(Iin, Ny, Nx, stepy, stepx, m, dm, fx, fy);
}	

// PARSEFLAG fft_filter FFT_Filter
void FFT_Filter(char *in)
{
	char *arg;
	char *opt;
	char *word;
	int Nx, Ny, stepy, stepx, m, dm;
	double fx, fy;
	
	if (!Iin.I)
	{
		AddErr(ERRNOIIMAGE);
		return;
	}
	
	// Parse options:
	word=malloc(strlen(in)*sizeof(char));
	
	GetOption(in, "Nx", word);
	Nx=atoi(word);	
	
	GetOption(in, "Ny", word);
	Ny=atoi(word);
	
	GetOption(in, "stepx", word);
	stepx=atoi(word);
	
	GetOption(in, "stepy", word);
	stepy=atoi(word);
	
	GetOption(in, "m", word);
	m=atoi(word);
	
	GetOption(in, "dm", word);
	dm=atoi(word);
	
	GetOption(in, "fx", word);
	fx=atof(word);
	
	GetOption(in, "fy", word);
	fy=atof(word);
	
	if (Iout.I)
		FreeImage(&Iout);
	Iout=FFT_PolynomalFilter(Iin, Ny, Nx, stepy, stepx, m, dm, fx, fy);
}	

#define MAXLINELEN 4096
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
		ParseComm(line);
		fgets(line, MAXLINELEN-1, f);
	}
	free(line);
	if (Iout.I)
		FreeImage(&Iout);
	if (Iin.I)
		FreeImage(&Iin);
	
}

void ParseCommandline(int argc, char **argv)
{
	int i=0;
	for (i=0;i<argc;i++)
		ParseComm(argv[i]);
	if (Iout.I)
		FreeImage(&Iout);
	if (Iin.I)
		FreeImage(&Iin);
}
