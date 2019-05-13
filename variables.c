#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "floatimage.h"
#include "filter.h"
#include "error.h"
typedef enum {FILTERVAR, FILTERSETVAR, IMAGEVAR} VarType;
typedef struct {
	filter F;
	filterset FS;
	image I;
	char *name;
	VarType T;
} Var;

#define BLOCK 10
Var *variables=NULL;
int Nvar=0;

void InitVars()
{
	if (variables)
		return; /* will only initialize once */
	variables=calloc(BLOCK, sizeof(Var));
	Nvar=0;
}
void ClearVars()
{
	int i;
	if (!variables)
		return;
	for (i=0;i<Nvar;i++)
	{
		switch(variables[i].T)
		{
			case FILTERVAR:
				FreeFilter(&(variables[i].F));
				break;
			case FILTERSETVAR:
				FreeFilterSet(&(variables[i].FS));
				break;
			case IMAGEVAR:
				FreeImage(&(variables[i].I));
				break;
			default:
				break;
		}
		free(variables[i].name);
	}
	free(variables);
	variables=NULL;
	Nvar=0;
}
void ListVars()
{
	int i=0;
	for (i=0;i<Nvar;i++)
	{
		printf("Variable \"%s\": ", variables[i].name);
		switch(variables[i].T)
		{
			case FILTERVAR:
				printf("Filter");
				break;
			case FILTERSETVAR:
				printf("Filter Set");
				break;
			case IMAGEVAR:
				printf("Image");
				break;
			default:
				break;
		}
		printf("\n");
	}
}


int lookupvar(char *name)
{
	int i=0;
	int lk;
	lk=strlen(name);
	for (i=0;i<Nvar;i++)
	{
		if (strlen(variables[i].name)==lk)
			if (strncmp(variables[i].name, name, lk)==0)
				return i;
	}
	return -1;	
}

int LookupFilter(char *name, filter *F)
{
	int i;
	if ((i=lookupvar(name))<0)
		return 0;
	else
	{
		if (variables[i].T==FILTERVAR)
		{
			*F=variables[i].F;
			return 1;
		}
	}
	return 0;
} 

int LookupFilterSet(char *name, filterset *FS)
{
	int i;
	if ((i=lookupvar(name))<0)
		return 0;
	else
	{
		if (variables[i].T==FILTERSETVAR)
		{
			*FS=variables[i].FS;
			return 1;
		}
	}
	return 0;
} 

int LookupImage(char *name, image *I)
{
	int i;
	if ((i=lookupvar(name))<0)
		return 0;
	else
	{
		if (variables[i].T==IMAGEVAR)
		{
			*I=variables[i].I;
			return 1;
		}
	}
	return 0;
} 

void AddFilter(char *name, filter F)
{
	int i;
	if ((i=lookupvar(name))<0)
	{
		variables[Nvar].name=name;
		variables[Nvar].F=F;
		variables[Nvar].T=FILTERVAR;
		Nvar++;
		if (Nvar%BLOCK==0)
			variables=realloc(variables, Nvar+BLOCK);
		return;
	}
	if (variables[i].T==FILTERVAR)	
	{
		FreeFilter(&(variables[i].F));
		variables[i].F=F;
	}
	else
	{
		// ERRORFLAG ERRVARWRONGTYPE  "variable exists as a different type"
		AddErr(ERRVARWRONGTYPE);
	}
}

void AddFilterSet(char *name, filterset FS)
{
	int i;
	if ((i=lookupvar(name))<0)
	{
		variables[Nvar].name=name;
		variables[Nvar].FS=FS;
		variables[Nvar].T=FILTERSETVAR;
		Nvar++;
		if (Nvar%BLOCK==0)
			variables=realloc(variables, Nvar+BLOCK);
		return;
	}
	if (variables[i].T==FILTERSETVAR)	
	{
		FreeFilterSet(&(variables[i].FS));
		variables[i].FS=FS;
	}
	else
		AddErr(ERRVARWRONGTYPE);
}

void AddImage(char *name, image I)
{
	int i;
	if ((i=lookupvar(name))<0)
	{
		variables[Nvar].name=name;
		variables[Nvar].I=I;
		variables[Nvar].T=IMAGEVAR;
		Nvar++;
		if (Nvar%BLOCK==0)
			variables=realloc(variables, Nvar+BLOCK);
		return;
	}
	if (variables[i].T==IMAGEVAR)	
	{
		FreeImage(&(variables[i].I));
		variables[i].I=I;
	}
	else
		AddErr(ERRVARWRONGTYPE);
}

void RMVar(char *name)
{
	int i;
	if ((i=lookupvar(name))>-1)
	{
		switch(variables[i].T)
		{
			case FILTERVAR:
				FreeFilter(&(variables[i].F));
				break;
			case FILTERSETVAR:
				FreeFilterSet(&(variables[i].FS));
				break;
			case IMAGEVAR:
				FreeImage(&(variables[i].I));
				break;
			default:
				break;
		}
		free(variables[i].name);
		i++;
		for (;i<Nvar;i++)
			variables[i-1]=variables[i];
		Nvar--;
		if (Nvar%BLOCK==BLOCK-1)
			variables=realloc(variables, Nvar+1);
	}
}


