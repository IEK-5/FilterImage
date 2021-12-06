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
// PARSEFLAG txtread TXTReadFile "I=<image-variable> file=<filename>"
void TXTReadFile(char *in)
{
	image I;
	char *name; 
	char *file; 
	
	name=malloc((strlen(in)+1)*sizeof(char));
	file=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "I", name);
	if (GetArg(in, "file", file))
		I=FloatimageTXTRead(file);
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", name);
		AddImage(name, I);
	}	
	else
		free(name);
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

// PARSEFLAG testimage MakeTestImage "I=<image-variable> N=<rows> M=<columns> m=<order>"
void MakeTestImage(char *in)
{	
	image I;
	char *name, *word;
	int N=10, M=10;
	double p=2;
	word=malloc((strlen(in)+1)*sizeof(char));
	name=malloc((strlen(in)+1)*sizeof(char));
	
	GetArg(in, "I", name);
	if (GetOption(in, "N", word))
		N=atoi(word);
	if (GetOption(in, "M", word))
		M=atoi(word);
	if (GetOption(in, "m", word))
		p=atof(word);
	I=TestImage(N, M, p);
	
	if (!ERRORSTATE)
	{
		printf("Defining image \"%s\"\n", name);
		AddImage(name, I);
	}
	else
		free(name);
	
	free(word);
}
void dmStringParser(char *dmstr, int **dmx, int **dmy, double **f, int *N)
{
	char *p, *q, *r;
	int j=1;
	
	/* string has the form <f1>x<n1>y<m1>+<f2>*x<n2>y<m2>+... */
	p=dmstr;
	while (*p)
	{
		if (*p=='+')
			j++;
		p++;
	}
	(*N)=0;
	(*dmx)=malloc(j*sizeof(int));
	(*dmy)=malloc(j*sizeof(int));
	(*f)=malloc(j*sizeof(double));
	
	p=dmstr;
	j=0;
	while (*p)
	{
		q=p;
		r=p;
		while (*q)
		{
			if (*q=='x')
			{
				*q='\0'; /* r points to a string <f>*/
				(*f)[j]=atof(r);
				q++;
				break;
			}
			q++;
		}
		if (!*q)
		{
			// ERRORFLAG ERRDERIVSTRINC  "incomplete filter spec"
			AddErr(ERRDERIVSTRINC);
			free(*dmx);
			free(*dmy);
			free(*f);
			*dmx=NULL;
			*dmy=NULL;
			*f=NULL;
			return;
		}
		r=q;
		while (*q)
		{
			if (*q=='y')
			{
				*q='\0'; /* r points to a string <n>*/
				(*dmx)[j]=atoi(r);
				q++;
				break;
			}
			q++;
		}
		if (!*q)
		{
			AddErr(ERRDERIVSTRINC);
			free(*dmx);
			free(*dmy);
			free(*f);
			*dmx=NULL;
			*dmy=NULL;
			*f=NULL;
			return;
		}
		r=q;
		while (*q)
		{
			if (*q=='+')
			{
				*q='\0'; /* r points to a string <m>*/
				(*dmy)[j]=atoi(r);
				q++;
				break;
			}
			q++;
		}
		if (!*q)
			(*dmy)[j]=atoi(r);
		p=q;
		j++;
		(*N)++;
	}
}
// PARSEFLAG makefilter MakeFilter "F=<filter-variable> nn=<int> ns=<int> ne=<int> nw=<int> m=<int> dm=<float>x<int>y<int>+... method=<[\'i\',\'q\',\'s\',\'p\']>"
void MakeFilter(char *in)
{
	filter F;	
	char *word;
	char *name;
	int nn=3, ns=3, nw=3, ne=3, m=2;
	char method='p';
	int *dmx=NULL, *dmy=NULL, N;
	double *f=NULL;
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
		
	if (GetOption(in, "method", word))
		method=word[0];	
	
	if (GetOption(in, "dm", word))
		dmStringParser(word, &dmx, &dmy, &f, &N);
	if (ERRORSTATE)
	{
		if (dmx)
			free(dmx);
		if (dmy)
			free(dmy);
		if (f)
			free(f);
		free(word);
		free(name);
		return;
	}
	
	if (ERRORSTATE)
	{
		free(word);
		free(name);
		return;
	}
	F=PartDeriv2DSum(nn, ns, nw, ne, m, dmx, dmy, f, N, method);
	if (!ERRORSTATE)
	{
		printf("Defining filter \"%s\"\n", name);
		AddFilter(name, F); /* do not free name! */
	}
	free(word);
	if (dmx)
		free(dmx);
	if (dmy)
		free(dmy);
	if (f)
		free(f);
}


// PARSEFLAG makefilterset MakeFilterSet "F=<filter-set-variable> nn=<int> ns=<int> ne=<int> nw=<int> m=<int> dm=<float>x<int>y<int>+... method=<[\'i\',\'q\',\'s\',\'p\']>"
void MakeFilterSet(char *in)
{
	filterset F;	
	char *word;
	char *name;
	int nn=3, ns=3, nw=3, ne=3, m=2;
	char method='p';
	int *dmx=NULL, *dmy=NULL, N;
	double *f=NULL;
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
	
	if (GetOption(in, "method", word))
		method=word[0];	
	
	if (GetOption(in, "dm", word))
		dmStringParser(word, &dmx, &dmy, &f, &N);
	
	if (ERRORSTATE)
	{
		if (dmx)
			free(dmx);
		if (dmy)
			free(dmy);
		if (f)
			free(f);
		free(word);
		free(name);
		return;
	}
	F=DerivOperatorSet2D(nn, ns, nw, ne, m, dmx, dmy, f, N,method);
	if (!ERRORSTATE)
	{
		printf("Defining filter set \"%s\"\n", name);
		AddFilterSet(name, F); /* do not free name! */
	}
	free(word);
	if (dmx)
		free(dmx);
	if (dmy)
		free(dmy);
	if (f)
		free(f);
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

// PARSEFLAG addimages AddIm "addimage I1=<image-variable> f1=<float> I2=<image-variable> f2=<float> I3=<image-variable>"
void AddIm(char *in)
{
	image I1, I2, I3;
	double f1, f2;
	char *name;
	name=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "I1", name);
	if (!LookupImage(name, &I1))
	{
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	GetArg(in, "I2", name);
	if (!LookupImage(name, &I2))
	{
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	GetArg(in, "f1", name);
	f1=atof(name);
	GetArg(in, "f2", name);
	f2=atof(name);
	GetArg(in, "I3", name);
	if (!ERRORSTATE)
	{
		I3=DupImage(I1);
		AddImages(I1, I2, f1, f2, I3);
		
		if (!ERRORSTATE)
		{
			printf("Defining image \"%s\"\n", name);
			AddImage(name, I3);
		}
		else
			FreeImage(&I3);
	}
	free(name);
}
// PARSEFLAG absimage AbsIm "absimage I=<image-variable>"
void AbsIm(char *in)
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
	if (!ERRORSTATE)
		AbsImage(I);
	free(name);
}

// PARSEFLAG clipimage ClipIm "clipimage I=<image-variable> th=<float> p=<int>"
void ClipIm(char *in)
{
	image I;
	int p;
	double th;
	char *name;
	name=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "I", name);
	if (!LookupImage(name, &I))
	{
		AddErr(ERRNOOIMAGE);
		free(name);
		return;
	}
	GetArg(in, "th", name);
	th=atof(name);
	GetArg(in, "p", name);
	p=atoi(name);
	if (!ERRORSTATE)
		ImageClip(I,th, p);
	free(name);
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


char * keyword_generator(const char *text, int state)
{
    static int list_index, len;
    char *name;

    if (!state) {
        list_index = 0;
        len = strlen(text);
    }

    while ((name = KeyTable[list_index++].key)) {
        if (strncmp(name, text, len) == 0) {
            return strdup(name);
        }
    }

    return NULL;
}


