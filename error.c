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
#include "error.h"
#include "errormessages.h" /* generated by gen_errorflags.sh, defines error messages for eachg error flag */

int ERRORSTATE=0;			/* single boolean indicating whether there were any errors */
int EXITSTATE=0;
char ERRORS[NERR] = {0};	/* zero initialized error flags */

void AddErr(int ERRFLAG)
{
	ERRORSTATE=1;
	EXITSTATE=1;
	if (ERRFLAG>NERR)
	{
		/* the recursive error :) */
		// ERRORFLAG ERROUTRANGE  "Error flag out of range!"
		AddErr(ERROUTRANGE);
	}
	else
		ERRORS[ERRFLAG]++;
}
int QueryErr(int ERRFLAG)
{
	if (ERRFLAG>NERR)
		AddErr(ERROUTRANGE);
	if (ERRORS[ERRFLAG])
		return 1;
	return 0;
}

void E_Messages()
{
	int i;
	for (i=0;i<NERR;i++)
		if (ERRORS[i])
			fprintf(stderr,"ERROR: %s (%dx)\n",  EMessages[i], ERRORS[i]);
}

void ResetErrors()
{
	int i;
	for (i=0;i<NERR;i++)
		ERRORS[i]=0;
	ERRORSTATE=0;
}
