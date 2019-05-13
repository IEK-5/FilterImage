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
		if (isspace(*end))
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
		if (!isspace(*end))
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
	int go=1;
	ParserFun fun;
	/* skip space chars */
	if (*in=='#')
		return;
	while (*in && go)
	{
		if (!isspace(*in))
		{
			go=0;
		}
		else
			in++;
	}
	if (!(*in))
		return;
	key=malloc((strlen(in)+1)*sizeof(char));
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

int GetOption(char *in, char *opt, char *word)
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
		*word='\0';
		free(opti);	
		return 0;
	}	
	GetWord(start+len-1, word);
	free(opti);	
	return 1;
}



/* parse routines
 * for each parse routine add a parseflag:
 * '// PARSEFLAG <keyword> <function-name>'
 * Note: you can create aliases by adding several lines like this with 
 * different keywords but the same function.
 * a parse routine must be of type void and take one character string as input
 * The gen_parseflags will take care of all the definitions in parsedef.h
 */
 
// PARSEFLAG imread ImRead
// PARSEFLAG imload ImRead
void ImRead(char *in)
{
	image I;
	char *name;
	name=malloc((strlen(in)+1)*sizeof(char));
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

// PARSEFLAG imsave ImSave
void ImSave(char *in)
{
	image I;
	char *name;
	name=malloc((strlen(in)+1)*sizeof(char));
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
	name=malloc((strlen(in)+1)*sizeof(char));
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
// PARSEFLAG fsetread FilterSetRead
// PARSEFLAG fsetload FilterSetRead
void FilterSetRead(char *in)
{
	filterset F;
	char *name;
	name=malloc((strlen(in)+1)*sizeof(char));
	in=GetWord(in, name); /* fetch the firt word of in, will be the name of the image */
	if (!in)
	{
		// ERRORFLAG ERRPREMEND  "premature end of input"
		AddErr(ERRPREMEND);
		free(name);
		return;
	}
	F=LoadFilterSet(in);
	if (!ERRORSTATE)
	{
		printf("Defining filterset \"%s\"\n", name);
		AddFilterSet(name, F);
	}		
}

// PARSEFLAG fsetsave FilterSetSave
void FilterSetSave(char *in)
{
	filterset F;
	char *name;
	name=malloc((strlen(in)+1)*sizeof(char));
	in=GetWord(in, name); /* fetch the firt word of in, will be the name of the image */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(name);
		return;
	}
	if (!LookupFilterSet(name, &F))
	{
		// ERRORFLAG ERRNOOIMAGE  "image not available"
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	free(name);
	if (!ERRORSTATE)
		SaveFilterSet(in, F);
}

// PARSEFLAG fread FilterRead
// PARSEFLAG fload FilterRead
void FilterRead(char *in)
{
	filter F;
	char *name;
	name=malloc((strlen(in)+1)*sizeof(char));
	in=GetWord(in, name); /* fetch the firt word of in, will be the name of the image */
	if (!in)
	{
		// ERRORFLAG ERRPREMEND  "premature end of input"
		AddErr(ERRPREMEND);
		free(name);
		return;
	}
	F=LoadFilter(in);
	if (!ERRORSTATE)
	{
		printf("Defining filter \"%s\"\n", name);
		AddFilter(name, F);
	}		
}

// PARSEFLAG fsave FilterSave
void FilterSave(char *in)
{
	filter F;
	char *name;
	name=malloc((strlen(in)+1)*sizeof(char));
	in=GetWord(in, name); /* fetch the firt word of in, will be the name of the image */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(name);
		return;
	}
	if (!LookupFilter(name, &F))
	{
		// ERRORFLAG ERRNOOIMAGE  "image not available"
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	free(name);
	if (!ERRORSTATE)
		SaveFilter(in, F);
}
// PARSEFLAG makefilter MakeFilter
void MakeFilter(char *in)
{
	filter F;	
	char *word;
	char *name;
	int nn=3, ns=3, nw=3, ne=3, dm=0, m=2;
	char method='p';
	double fx=1, fy=1;
	word=malloc((strlen(in)+1)*sizeof(char));
	name=malloc((strlen(in)+1)*sizeof(char));
	in=GetWord(in, name); /* fetch the first word of in, will be the name of the input image */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(word);
		return;
	}
	
	if (GetOption(in, "nn", word))
		nn=atoi(word);	
	
	if (GetOption(in, "ns", word))
		ns=atoi(word);
	
	if (GetOption(in, "ne", word))
		ne=atoi(word);
	
	if (GetOption(in, "nw", word))
		nw=atoi(word);
	
	if (GetOption(in, "m", word))
		m=atoi(word);
	
	if (GetOption(in, "dm", word))
		dm=atoi(word);
	
	if (GetOption(in, "fx", word))
		fx=atof(word);
	
	if (GetOption(in, "fy", word))
		fy=atof(word);
		
	if (GetOption(in, "method", word))
		method=word[0];	
	
	
	if (ERRORSTATE)
	{
		free(word);
		free(name);
		return;
	}
	F=PartDeriv2D(nn, ns, nw, ne, dm, m, fx, fy, method);
	if (!ERRORSTATE)
	{
		printf("Defining filter \"%s\"\n", name);
		AddFilter(name, F); /* do not free name! */
	}
	free(word);
}

// PARSEFLAG makefilterset MakeFilterSet
void MakeFilterSet(char *in)
{
	filterset F;	
	char *word;
	char *name;
	int nn=3, ns=3, nw=3, ne=3, dm=0, m=2;
	char method='p';
	double fx=1, fy=1;
	word=malloc((strlen(in)+1)*sizeof(char));
	name=malloc((strlen(in)+1)*sizeof(char));
	in=GetWord(in, name); /* fetch the first word of in, will be the name of the input image */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(word);
		return;
	}
	
	if (GetOption(in, "nn", word))
		nn=atoi(word);	
	
	if (GetOption(in, "ns", word))
		ns=atoi(word);
	
	if (GetOption(in, "ne", word))
		ne=atoi(word);
	
	if (GetOption(in, "nw", word))
		nw=atoi(word);
	
	if (GetOption(in, "m", word))
		m=atoi(word);
	
	if (GetOption(in, "dm", word))
		dm=atoi(word);
	
	if (GetOption(in, "fx", word))
		fx=atof(word);
	
	if (GetOption(in, "fy", word))
		fy=atof(word);
		
	if (GetOption(in, "method", word))
		method=word[0];	
	
	
	if (ERRORSTATE)
	{
		free(word);
		free(name);
		return;
	}
	F=DerivOperatorSet2D(nn, ns, nw, ne, dm, m, fx, fy,method);
	if (!ERRORSTATE)
	{
		printf("Defining filter set \"%s\"\n", name);
		AddFilterSet(name, F); /* do not free name! */
	}
	free(word);
}


// PARSEFLAG plain_applyfilterset PlainApplyFilterSet
void PlainApplyFilterSet(char *in)
{
	char *word;
	char *name;
	int stepy=1, stepx=1;
	image I, Iout;
	filterset F;
	
	word=malloc((strlen(in)+1)*sizeof(char));
	name=malloc((strlen(in)+1)*sizeof(char));
	in=GetWord(in, word); /* fetch the first word of in, will be the name of the input image */
	in=GetWord(in, name); /* fetch the second word of in, will be the name of the output image */
	if (!LookupImage(word, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(word);
		return;
	}
	in=GetWord(in, word); /* fetch the filter */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(word);
		return;
	}
	if (!LookupFilterSet(word, &F))
	{
		// ERRORFLAG ERRNOFILTERSET  "filterset not available"
		AddErr(ERRNOFILTERSET);
		free(word);
		return;
	}
	// Parse options:
	
	if (GetOption(in, "stepx", word))
		stepx=atoi(word);
	
	if (GetOption(in, "stepy", word))
		stepy=atoi(word);
	
	if (ERRORSTATE)
	{
		free(word);
		free(name);
		return;
	}
	Iout=PL_ApplyFilter(I, stepy, stepx,  F);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", name);
		AddImage(name, Iout);
	}
	free(word);
}	

// PARSEFLAG applyfilterset ApplyFilterSet
void ApplyFilterSet(char *in)
{
	char *word;
	char *name;
	int stepy=1, stepx=1;
	image I, Iout;
	filterset F;
	
	word=malloc((strlen(in)+1)*sizeof(char));
	name=malloc((strlen(in)+1)*sizeof(char));
	in=GetWord(in, word); /* fetch the first word of in, will be the name of the input image */
	in=GetWord(in, name); /* fetch the second word of in, will be the name of the output image */
	if (!LookupImage(word, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(word);
		return;
	}
	in=GetWord(in, word); /* fetch the filter */
	if (!in)
	{
		AddErr(ERRPREMEND);
		free(word);
		return;
	}
	if (!LookupFilterSet(word, &F))
	{
		// ERRORFLAG ERRNOFILTERSET  "filterset not available"
		AddErr(ERRNOFILTERSET);
		free(word);
		return;
	}
	// Parse options:
	
	if (GetOption(in, "stepx", word))
		stepx=atoi(word);
	
	if (GetOption(in, "stepy", word))
		stepy=atoi(word);
	
	if (ERRORSTATE)
	{
		free(word);
		free(name);
		return;
	}
	Iout=ApplyFilter(I, stepy, stepx,  F);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", name);
		AddImage(name, Iout);
	}
	free(word);
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
	name=malloc((strlen(in)+1)*sizeof(char));
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
	int Nx=3, Ny=3, stepy=1, stepx=1, m=2, dm=0;
	double fx=1, fy=1;
	char method='p';
	image I, Iout;
	word=malloc((strlen(in)+1)*sizeof(char));
	name=malloc((strlen(in)+1)*sizeof(char));
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
	
	if (GetOption(in, "Nx", word))
		Nx=atoi(word);	
	
	if (GetOption(in, "Ny", word))
		Ny=atoi(word);
	
	if (GetOption(in, "stepx", word))
		stepx=atoi(word);
	
	if (GetOption(in, "stepy", word))
		stepy=atoi(word);
	
	if (GetOption(in, "m", word))
		m=atoi(word);
	
	if (GetOption(in, "dm", word))
		dm=atoi(word);
	
	if (GetOption(in, "fx", word))
		fx=atof(word);
	
	if (GetOption(in, "fy", word))
		fy=atof(word);
	
	if (GetOption(in, "method", word))
		method=word[0];	
	if (ERRORSTATE)
	{
		free(word);
		free(name);
		return;
	}
	Iout=PL_PolynomalFilter(I, Ny, Nx, stepy, stepx, m, dm, fx, fy,method);
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
	int Nx=3, Ny=3, stepy=1, stepx=1, m=2, dm=0;
	double fx=1, fy=1;
	char method='p';
	
	image I, Iout;
	word=malloc((strlen(in)+1)*sizeof(char));
	name=malloc((strlen(in)+1)*sizeof(char));
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
	
	if (GetOption(in, "Nx", word))
		Nx=atoi(word);	
	
	if (GetOption(in, "Ny", word))
		Ny=atoi(word);
	
	if (GetOption(in, "stepx", word))
		stepx=atoi(word);
	
	if (GetOption(in, "stepy", word))
		stepy=atoi(word);
	
	if (GetOption(in, "m", word))
		m=atoi(word);
	
	if (GetOption(in, "dm", word))
		dm=atoi(word);
	
	if (GetOption(in, "fx", word))
		fx=atof(word);
	
	if (GetOption(in, "fy", word))
		fy=atof(word);
	
	if (GetOption(in, "method", word))
		method=word[0];	
	if (ERRORSTATE)
	{
		free(word);
		free(name);
		return;
	}
	Iout=PolynomalFilter(I, Ny, Nx, stepy, stepx, m, dm, fx, fy,method);
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
	int Nx=3, Ny=3, stepy=1, stepx=1, m=2, dm=0;
	double fx=1, fy=1;
	char method='p';
	
	image I, Iout;
	word=malloc((strlen(in)+1)*sizeof(char));
	name=malloc((strlen(in)+1)*sizeof(char));
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
	if (GetOption(in, "Nx", word))
		Nx=atoi(word);	
	
	if (GetOption(in, "Ny", word))
		Ny=atoi(word);
	
	if (GetOption(in, "stepx", word))
		stepx=atoi(word);
	
	if (GetOption(in, "stepy", word))
		stepy=atoi(word);
	
	if (GetOption(in, "m", word))
		m=atoi(word);
	
	if (GetOption(in, "dm", word))
		dm=atoi(word);
	
	if (GetOption(in, "fx", word))
		fx=atof(word);
	
	if (GetOption(in, "fy", word))
		fy=atof(word);
		
	if (GetOption(in, "method", word))
		method=word[0];	
	
	if (ERRORSTATE)
	{
		free(word);
		free(name);
		return;
	}
	Iout=FFT_PolynomalFilter(I, Ny, Nx, stepy, stepx, m, dm, fx, fy,method);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", name);
		AddImage(name, Iout);
	}
	free(word);
}	


