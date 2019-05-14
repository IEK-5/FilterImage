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
	int squoted=0, dquoted=0;
	int go=1;
	while (*end && go)
	{
		if (*end=='\'')
			squoted=!squoted;
		if (*end=='\"')
			dquoted=!dquoted;
			
		if (!(squoted||dquoted)&&(isspace(*end)))
		{
			go=0;
		}
		else 
		{
			if ((*end!='\'')||(*end=='\"'))
			{
				*out=*end;
				out++;
			}
		}
		end++;
	}
	*out='\0';
	if (squoted||dquoted)
	{
		// ERRORFLAG ERRUNMQUO  "unmatched quotes"
		AddErr(ERRUNMQUO);
	}
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


int ParseComm(char *in)
{
	char *key;
	char *arg;
	int go=1;
	ParserFun fun;
	/* skip space chars */
	if (*in=='#')
		return 0;
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
		return 0;
	key=malloc((strlen(in)+1)*sizeof(char));
	arg=GetWord(in, key);
	if (strncmp(key, "exit", 5)==0)
	{
		free(key);
		return 1;
	}
		
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
	return 0;
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

int GetArg(char *in, char *opt, char *word)
{
	if (!GetOption(in, opt, word))
	{		
		// ERRORFLAG ERRMISSINGARG  "Mandatory argument missing"
		AddErr(ERRMISSINGARG);
		return 0;
	}
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
 
// PARSEFLAG imread ImRead "I=<image-variable> file=<filename>"
// PARSEFLAG imload ImRead "I=<image-variable> file=<filename>"
void ImRead(char *in)
{
	image I;
	char *name; 
	char *file; 
	
	name=malloc((strlen(in)+1)*sizeof(char));
	file=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "I", name);
	if (GetArg(in, "file", file))
		I=FloatimageRead(file);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", name);
		AddImage(name, I);
	}	
	else
		free(name);
	free(file);
}

// PARSEFLAG imsave ImSave "I=<image-variable> file=<filename>"
void ImSave(char *in)
{
	image I;
	char *name;
	char *file; 
	
	name=malloc((strlen(in)+1)*sizeof(char));
	file=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "I", name);
	GetArg(in, "file", file);
	if (!LookupImage(name, &I))
	{
		// ERRORFLAG ERRNOOIMAGE  "image not available"
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	free(name);
	if (!ERRORSTATE)
		FloatimageWrite(file, I, 1,0,1);
	free(file);
}

// PARSEFLAG show Show "I=<image-variable>"
void Show(char *in)
{
	image I;
	char *name;
	
	name=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "I", name);
	if (!LookupImage(name, &I))
	{
		// ERRORFLAG ERRNOOIMAGE  "image not available"
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	free(name);
	if (!ERRORSTATE)
		FloatimageDisplay(I, 1, 0, 1);
}

// PARSEFLAG txtsave TXTSaveFile "I=<image-variable> file=<filename>"
void TXTSaveFile(char *in)
{
	image I;
	char *name;
	char *file; 
	
	name=malloc((strlen(in)+1)*sizeof(char));
	file=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "I", name);
	GetArg(in, "file", file);
	if (!LookupImage(name, &I))
	{
		// ERRORFLAG ERRNOOIMAGE  "image not available"
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	free(name);
	if (!ERRORSTATE)
	{
		if (strncmp(file, "-", 2)==0)
			Floatimage2stdout(I);
		else
			FloatimageTXTWrite(file, I);
	}
	free(file);
}
// PARSEFLAG fsetread FilterSetRead "F=<filter-set-variable> file=<filename>"
// PARSEFLAG fsetload FilterSetRead "F=<filter-set-variable> file=<filename>"
void FilterSetRead(char *in)
{
	filterset F;
	char *name;
	char *file; 
	
	name=malloc((strlen(in)+1)*sizeof(char));
	file=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "F", name);
	GetArg(in, "file", file);
	
	if (!ERRORSTATE)
		F=LoadFilterSet(file);
	if (!ERRORSTATE)
	{
		printf("Defining filterset \"%s\"\n", name);
		AddFilterSet(name, F);
	}
	else
		free(name);
	free(file);
}

// PARSEFLAG fsetsave FilterSetSave "F=<filter-set-variable> file=<filename>"
void FilterSetSave(char *in)
{
	filterset F;
	char *name;
	char *file; 
	
	name=malloc((strlen(in)+1)*sizeof(char));
	file=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "F", name);
	GetArg(in, "file", file);
	if (!LookupFilterSet(name, &F))
	{
		// ERRORFLAG ERRNOFSET  "filter set not available"
		AddErr(ERRNOFSET);
		free(name);
		return;
	}
	free(name);
	if (!ERRORSTATE)
		SaveFilterSet(file, F);
	free(file);
}

// PARSEFLAG fread FilterRead "F=<filter-variable> file=<filename>"
// PARSEFLAG fload FilterRead "F=<filter-variable> file=<filename>"
void FilterRead(char *in)
{
	filter F;
	char *name;
	char *file; 
	
	name=malloc((strlen(in)+1)*sizeof(char));
	file=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "F", name);
	GetArg(in, "file", file);
	
	if (!ERRORSTATE)
		F=LoadFilter(file);
	if (!ERRORSTATE)
	{
		printf("Defining filter \"%s\"\n", name);
		AddFilter(name, F);
	}	
	else
		free(name);
	free(file);
}

// PARSEFLAG fsave FilterSave "F=<filter-variable> file=<filename>"
void FilterSave(char *in)
{
	filter F;
	char *name;
	char *file; 
	
	name=malloc((strlen(in)+1)*sizeof(char));
	file=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "F", name);
	GetArg(in, "file", file);
	
	if (!LookupFilter(name, &F))
	{
		// ERRORFLAG ERRNOOIMAGE  "image not available"
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	free(name);
	if (!ERRORSTATE)
		SaveFilter(file, F);
	free(file);
}
// PARSEFLAG makefilter MakeFilter "F=<filter-variable> nn=<int> ns=<int> ne=<int> nw=<int> m=<int> dm=<int> fx=<float> fy=<float> method=<char>"
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
	GetArg(in, "F", name);
	
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

// PARSEFLAG makefilterset MakeFilterSet "F=<filter-set-variable> nn=<int> ns=<int> ne=<int> nw=<int> m=<int> dm=<int> fx=<float> fy=<float> method=<char>"
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
	GetArg(in, "F", name);
	
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


// PARSEFLAG plain_applyfilterset PlainApplyFilterSet "Iin=<input-image-variable> Iout=<output-image-variable> F=<filter-set-variable> stepx=<int> stepy=<int>"
void PlainApplyFilterSet(char *in)
{
	char *word;
	int stepy=1, stepx=1;
	image I, Iout;
	filterset F;
	char *iname, *oname, *fname;
	
	word=malloc((strlen(in)+1)*sizeof(char));
	iname=malloc((strlen(in)+1)*sizeof(char));
	oname=malloc((strlen(in)+1)*sizeof(char));
	fname=malloc((strlen(in)+1)*sizeof(char));
	
	GetArg(in, "Iin", iname);
	GetArg(in, "Iout", oname);
	GetArg(in, "F", fname);
	
	if (!LookupImage(iname, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(iname);
		free(oname);
		free(fname);
		free(word);
		return;
	}
	
	if (!LookupFilterSet(fname, &F))
	{
		// ERRORFLAG ERRNOFILTERSET  "filterset not available"
		AddErr(ERRNOFILTERSET);
		free(iname);
		free(oname);
		free(fname);
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
		free(iname);
		free(oname);
		free(fname);
		free(word);
		return;
	}
	Iout=PL_ApplyFilter(I, stepy, stepx,  F);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", oname);
		AddImage(oname, Iout);
	}
	else
		free(oname);
	free(iname);
	free(fname);
	free(word);
}	

// PARSEFLAG applyfilterset ApplyFilterSet "Iin=<input-image-variable> Iout=<output-image-variable> F=<filter-set-variable> stepx=<int> stepy=<int>"
void ApplyFilterSet(char *in)
{
	char *word;
	int stepy=1, stepx=1;
	image I, Iout;
	filterset F;
	
	char *iname, *oname, *fname;
	
	word=malloc((strlen(in)+1)*sizeof(char));
	iname=malloc((strlen(in)+1)*sizeof(char));
	oname=malloc((strlen(in)+1)*sizeof(char));
	fname=malloc((strlen(in)+1)*sizeof(char));
	
	GetArg(in, "Iin", iname);
	GetArg(in, "Iout", oname);
	GetArg(in, "F", fname);
	
	if (!LookupImage(iname, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(iname);
		free(oname);
		free(fname);
		free(word);
		return;
	}
	
	if (!LookupFilterSet(fname, &F))
	{
		// ERRORFLAG ERRNOFILTERSET  "filterset not available"
		AddErr(ERRNOFILTERSET);
		free(iname);
		free(oname);
		free(fname);
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
		free(iname);
		free(oname);
		free(fname);
		free(word);
		return;
	}
	Iout=ApplyFilter(I, stepy, stepx,  F);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", oname);
		AddImage(oname, Iout);
	}
	else
		free(oname);
	free(iname);
	free(fname);
	free(word);
}	

// PARSEFLAG fft_applyfilter FFTApplyFilter "Iin=<input-image-variable> Iout=<output-image-variable> F=<filter-variable> stepx=<int> stepy=<int>"
void FFTApplyFilter(char *in)
{
	char *word;
	int stepy=1, stepx=1;
	image I, Iout;
	filter F;
	
	char *iname, *oname, *fname;
	
	word=malloc((strlen(in)+1)*sizeof(char));
	iname=malloc((strlen(in)+1)*sizeof(char));
	oname=malloc((strlen(in)+1)*sizeof(char));
	fname=malloc((strlen(in)+1)*sizeof(char));
	
	GetArg(in, "Iin", iname);
	GetArg(in, "Iout", oname);
	GetArg(in, "F", fname);
	
	if (!LookupImage(iname, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(iname);
		free(oname);
		free(fname);
		free(word);
		return;
	}
	
	if (!LookupFilter(fname, &F))
	{
		// ERRORFLAG ERRNOFILTERSET  "filterset not available"
		AddErr(ERRNOFILTERSET);
		free(iname);
		free(oname);
		free(fname);
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
		free(iname);
		free(oname);
		free(fname);
		free(word);
		return;
	}
	Iout=FFT_ApplyFilter(I, stepy, stepx,  F);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", oname);
		AddImage(oname, Iout);
	}
	else
		free(oname);
	free(iname);
	free(fname);
	free(word);
}	

// PARSEFLAG who Who "who"
void Who(char *in)
{
	ListVars();
}
// PARSEFLAG ELtoVj ELtoVj "ELtoVj I=<image-variable> T=<float>"
void ELtoVj(char *in)
{
	image I;
	char *name;
	name=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "I", name);
	if (!LookupImage(name, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	GetArg(in, "T", name);
	if (!ERRORSTATE)	
		EL2Vj(I, atof(name)); // convert to junction voltages
	free(name);
}


// PARSEFLAG plain_filter PlainFilter "Iin=<input-image-variable> Iout=<output-image-variable> Nx=<int> Ny=<int> stepx=<int> stepy=<int> m=<int> dm=<int> fx=<float> fy=<float> method=<char>"
void PlainFilter(char *in)
{
	char *word;
	int Nx=3, Ny=3, stepy=1, stepx=1, m=2, dm=0;
	double fx=1, fy=1;
	char method='p';
	image I, Iout;
	char *iname, *oname;
	word=malloc((strlen(in)+1)*sizeof(char));
	iname=malloc((strlen(in)+1)*sizeof(char));
	oname=malloc((strlen(in)+1)*sizeof(char));
	
	GetArg(in, "Iin", iname);
	GetArg(in, "Iout", oname);
	if (!LookupImage(iname, &I))
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
		free(iname);
		free(oname);
		return;
	}
	Iout=PL_PolynomalFilter(I, Ny, Nx, stepy, stepx, m, dm, fx, fy,method);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", oname);
		AddImage(oname, Iout);
	}
	else
		free(oname);
	free(iname);
	free(word);
}	
	
// PARSEFLAG filter Filter "Iin=<input-image-variable> Iout=<output-image-variable> Nx=<int> Ny=<int> stepx=<int> stepy=<int> m=<int> dm=<int> fx=<float> fy=<float> method=<char>"
void Filter(char *in)
{
	char *word;
	int Nx=3, Ny=3, stepy=1, stepx=1, m=2, dm=0;
	double fx=1, fy=1;
	char method='p';	
	image I, Iout;
	char *iname, *oname;
	word=malloc((strlen(in)+1)*sizeof(char));
	iname=malloc((strlen(in)+1)*sizeof(char));
	oname=malloc((strlen(in)+1)*sizeof(char));
	
	GetArg(in, "Iin", iname);
	GetArg(in, "Iout", oname);
	if (!LookupImage(iname, &I))
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
		free(iname);
		free(oname);
		return;
	}	
	Iout=PolynomalFilter(I, Ny, Nx, stepy, stepx, m, dm, fx, fy,method);
	if (!ERRORSTATE)
	{		
		printf("Defining image \"%s\"\n", oname);
		AddImage(oname, Iout);
	}
	else
		free(oname);
	free(iname);
	free(word);
}	

// PARSEFLAG fft_filter FFT_Filter "Iin=<input-image-variable> Iout=<output-image-variable> Nx=<int> Ny=<int> stepx=<int> stepy=<int> m=<int> dm=<int> fx=<float> fy=<float> method=<char>"
void FFT_Filter(char *in)
{
	char *word;
	int Nx=3, Ny=3, stepy=1, stepx=1, m=2, dm=0;
	double fx=1, fy=1;
	char method='p';	
	image I, Iout;
	char *iname, *oname;
	word=malloc((strlen(in)+1)*sizeof(char));
	iname=malloc((strlen(in)+1)*sizeof(char));
	oname=malloc((strlen(in)+1)*sizeof(char));
	
	GetArg(in, "Iin", iname);
	GetArg(in, "Iout", oname);
	if (!LookupImage(iname, &I))
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
		free(iname);
		free(oname);
		return;
	}	
	Iout=FFT_PolynomalFilter(I, Ny, Nx, stepy, stepx, m, dm, fx, fy,method);
	if (!ERRORSTATE)
	{		
		printf("Defining image \"%s\"\n", oname);
		AddImage(oname, Iout);
	}
	else
		free(oname);
	free(iname);
	free(word);
}

void SplitWords(char *in, char *ident)
{
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	while(*in)
	{
		in=GetWord(in, word);
		printf("%s%s\n", ident, word);
	}
	free(word);
}

// PARSEFLAG help Help "[-l/command]"
void Help(char *in)
{
	int i=0;
	int lk;
	lk=strlen(in);
	if ((strncmp(in,"-l", 3)==0)||(!*in))
	{
		while(Usage[i])
		{
			printf("Command %s:\n",KeyTable[i].key);
			SplitWords(Usage[i], "\t");
			printf("\n");
			i++;
		}
		return;
	}
			
	while (KeyTable[i].key)
	{
		if (strlen(KeyTable[i].key)==lk)
			if (strncmp(KeyTable[i].key, in, lk)==0)
			{
				printf("Command %s:\n",KeyTable[i].key);
				SplitWords(Usage[i], "\t");
				printf("\n");
				return;
			}
		i++;
	}
	printf("Command %s not known\n",in);
}



