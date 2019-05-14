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
FILE="parsedef.h"
echo "/* File generated by gen_parserflags.sh, do not edit by hand */" >  $FILE
NFLAGS=0;

# first collect all error flags in one file
for s in $@
do
	echo Collecting error flags from $s
	egrep -o '^// PARSEFLAG.*' $s | sort | uniq >>tmpflags
done

# create the parse flag defines
awk '{print "void "$4"(char *in);"}' tmpflags | sort | uniq>>$FILE

echo "typedef void (*ParserFun)(char *in);">>$FILE
echo "typedef struct {">>$FILE
echo "	char *key;">>$FILE
echo "	ParserFun fun;">>$FILE
echo "} KeyWord;">>$FILE

echo "const KeyWord KeyTable[] = {">>$FILE
awk '{print "\t{\""$3"\", &"$4"},"}' tmpflags>>$FILE
echo "	{NULL, NULL}">>$FILE
echo "};">>$FILE


echo "char *Usage[] = {"  >>  $FILE
egrep -o '// PARSEFLAG.*' tmpflags|egrep -o '\".*\"' |awk '{print "	"$0","}' >>  $FILE
echo "	NULL"  >>  $FILE
echo "};"   >>  $FILE


rm tmpflags
