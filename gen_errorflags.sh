#!/bin/bash

# Generate error flags and messages
# to be used in conjunction with error.c and error.h
# searches the sources for commentlines following the pattern:
# // ERRORFLAG <FLAGNAME> "<message string>"
#
# it generates a macro by the name FLAGNAME and adds a string to the error string array
# the defines are in errorflags.h"
# the message array in "errormessages.h"
# the errorflags.h should be included in all sources that throw errors
# the  "errormessages.h" is a private include for the error.c

echo "/* File generated by gen_errorflags.sh, do not edit by hand */" >  "errorflags.h"
NERR=0;

# first collect all error flags in one file
for s in $@
do
	echo Collecting error flags from $s
	egrep -o '// ERRORFLAG.*' $s >>tmperrflags
done

# create the error flag defines
awk '{print "#define",$3,NR-1}' tmperrflags>>"errorflags.h"
NERR=$(wc -l tmperrflags | awk '{print $1}')
	
	
echo "#define NERR $NERR" >>"errorflags.h"


echo "/* File generated by gen_errorflags.sh, do not edit by hand */" >  "errormessages.h"
echo "char *EMessages[] = {"  >>"errormessages.h"
egrep -o '// ERRORFLAG.*' tmperrflags|egrep -o '\".*\"' |awk '{print "	"$0","}' >>"errormessages.h"

echo "	NULL"  >>"errormessages.h"
echo "};"  >>"errormessages.h"


rm tmperrflags
