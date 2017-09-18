#! /bin/bash

set -e 
set -o pipefail

if [ "$#" -ne 4 ]; then
 	echo "args are : \n- inputfile (to sort)\n- outputfile (sorted)\n- number of lines of the inputfile\n- number of reducers required"
	exit
fi

INPUTFILE=$1
OUTPUTFILE=$2
NB_LINES=$3
NB_REDUCERS=$4


hdfs dfs -rm -r $INPUTFILE
hdfs dfs -copyFromLocal $INPUTFILE
hdfs dfs -rm -r $OUTPUTFILE
yarn jar target/sort-mapreduce-0.0.1.jar $INPUTFILE  $OUTPUFILE $NB_LINES $NB_REDUCERS
hdfs dfs -copyToLocal $OUTPUTFILE tests/$OUTPUTFILE
