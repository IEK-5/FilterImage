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
#include "variables.h"

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
	go=1; /* skip blank chars */
	while (*end && go)
	{
		if (!isblank(*end))
		{
			go=0;
		}
		else
			end++;
	}
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
	image I;
	char *name;
	name=malloc(strlen(in)*sizeof(char));
	in=GetWord(in, name); /* fetch the firt word of in, will be the name of the image */
	if (!in)
	{
		// ERRORFLAG ERRPREMEND  "premature end of input"
		AddErr(ERRPREMEND);
		free(name);
		return;
	}
	I=FloatimageRead(in);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", name);
		AddImage(name, I);
	}		
}

// PARSEFLAG save SaveFile
void SaveFile(char *in)
{
	image I;
	char *name;
	name=malloc(strlen(in)*sizeof(char));
	in=GetWord(in, name); /* fetch the firt word of in, will be the name of the image */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(name);
		return;
	}
	if (!LookupImage(name, &I))
	{
		// ERRORFLAG ERRNOOIMAGE  "image not available"
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	free(name);
	if (!ERRORSTATE)
		FloatimageWrite(in, I, 1,0,1);
}

// PARSEFLAG txtsave TXTSaveFile
void TXTSaveFile(char *in)
{
	image I;
	char *name;
	name=malloc(strlen(in)*sizeof(char));
	in=GetWord(in, name); /* fetch the firt word of in, will be the name of the image */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(name);
		return;
	}
	if (!LookupImage(name, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	free(name);
	if (ERRORSTATE)
		return;
	if (strncmp(in, "-", 2)==0)
		Floatimage2stdout(I);
	else
		FloatimageTXTWrite(in, I);
}

// PARSEFLAG who Who
void Who(char *in)
{
	ListVars();
}
// PARSEFLAG ELtoVj ELtoVj
void ELtoVj(char *in)
{
	image I;
	char *name;
	name=malloc(strlen(in)*sizeof(char));
	in=GetWord(in, name); /* fetch the firt word of in, will be the name of the image */
	if (!LookupImage(name, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	
	EL2Vj(I, atof(name)); // convert to junction voltages
	free(name);
}

/* todo split in makefilter, makefilterset, and applyfilterset, plain_applyfilterset, fft_applyfilter */

// PARSEFLAG plain_filter PlainFilter
void PlainFilter(char *in)
{
	char *word;
	char *name;
	int Nx, Ny, stepy, stepx, m, dm;
	double fx, fy;	
	image I, Iout;
	word=malloc(strlen(in)*sizeof(char));
	name=malloc(strlen(in)*sizeof(char));
	in=GetWord(in, word); /* fetch the first word of in, will be the name of the input image */
	in=GetWord(in, name); /* fetch the second word of in, will be the name of the output image */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(word);
		return;
	}
	if (!LookupImage(word, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(word);
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
	if (ERRORSTATE)
		return;
	dm=atoi(word);
	
	GetOption(in, "fx", word);
	fx=atof(word);
	
	GetOption(in, "fy", word);
	fy=atof(word);
	
	if (ERRORSTATE)
	{
		free(word);
		return;
	}
	Iout=PL_PolynomalFilter(I, Ny, Nx, stepy, stepx, m, dm, fx, fy);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", name);
		AddImage(name, Iout);
	}
	free(word);
}	
	
// PARSEFLAG filter Filter
void Filter(char *in)
{
	char *word;
	char *name;
	int Nx, Ny, stepy, stepx, m, dm;
	double fx, fy;
	
	image I, Iout;
	word=malloc(strlen(in)*sizeof(char));
	name=malloc(strlen(in)*sizeof(char));
	in=GetWord(in, word); /* fetch the first word of in, will be the name of the input image */
	in=GetWord(in, name); /* fetch the second word of in, will be the name of the output image */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(word);
		return;
	}
	if (!LookupImage(word, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(word);
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
	
	if (ERRORSTATE)
	{
		free(word);
		return;
	}
	Iout=PolynomalFilter(I, Ny, Nx, stepy, stepx, m, dm, fx, fy);
	if (!ERRORSTATE)
	{		
		printf("Defining image \"%s\"\n", name);
		AddImage(name, Iout);
	}
	free(word);
}	

// PARSEFLAG fft_filter FFT_Filter
void FFT_Filter(char *in)
{
	char *word;
	char *name;
	int Nx, Ny, stepy, stepx, m, dm;
	double fx, fy;	
	image I, Iout;
	word=malloc(strlen(in)*sizeof(char));
	name=malloc(strlen(in)*sizeof(char));
	in=GetWord(in, word); /* fetch the first word of in, will be the name of the input image */
	in=GetWord(in, name); /* fetch the second word of in, will be the name of the output image */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(word);
		return;
	}
	if (!LookupImage(word, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(word);
		return;
	}
	
	// Parse options:	
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
	
	if (ERRORSTATE)
	{
		free(word);
		return;
	}
	Iout=FFT_PolynomalFilter(I, Ny, Nx, stepy, stepx, m, dm, fx, fy);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", name);
		AddImage(name, Iout);
	}
	free(word);
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
	InitVars();
	line=malloc(MAXLINELEN*sizeof(char));
    fgets(line, MAXLINELEN-1, f);
	while(feof(f)==0)
	{
		ParseComm(line);
		fgets(line, MAXLINELEN-1, f);
	}
	free(line);
	ClearVars();
}

void ParseCommandline(int argc, char **argv)
{
	int i=0;
	InitVars();
	for (i=0;i<argc;i++)
		ParseComm(argv[i]);
	ClearVars();
}
