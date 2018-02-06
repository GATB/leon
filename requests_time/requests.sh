#!/bin/bash

if [ "$#" -ne 5 ]; then
 	echo "args are : \n- requests file\n- input file\n- outputfile\n- stdout file\n- stderr file\n"
	exit 1
fi

#arguments
REQUESTS_FILE=$1
INPUT_FILE=$2
OUTPUT_FILE=$3
STDOUT=$4
STDERR=$5

../leon -file "$INPUT_FILE" -output "$OUTPUT_FILE" -c -r < "$REQUESTS_FILE" > "$STDOUT" 2> "$STDERR"

