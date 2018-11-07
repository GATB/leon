#! /bin/bash

set -e 
#set -o pipefail

ABSOLUTE_PATH_TO_REPOSITORY="/home/tbraquel/NGS/gatb-tool-leon"

if [ "$ABSOLUTE_PATH_TO_REPOSITORY" = "" ]; then
	echo "please give the absolute path to the leon repository in the script launch_mapred_sort.sh, in the variable ABSOLUTE_PATH_TO_REPOSITORY"
	exit 1
fi


if [ "$#" -ne 5 ]; then
 	echo "args are : \n- path to input/output files\n- inputfile (to sort)\n- outputfile (sorted)\n- number of lines of the inputfile\n- number of reducers required"
	exit 1
fi

#TODO
#create tmp directory for users
HDFS_PATH="hdfs://mycluster/user/tbraquel/"

PATH_TO_FILES=$1
INPUTFILE=$2
OUTPUTFILE=$3
NB_LINES=$4
NB_REDUCERS=$5

#TESTS
echo "test args : "
echo "arg 1 = " $PATH_TO_FILES
echo "arg 2 = " $INPUTFILE
echo "arg 3 = " $OUTPUTFILE
echo "arg 4 = " $NB_LINES
echo "arg 5 = " $NB_REDUCERS

hdfs dfs -rm -r -f $INPUTFILE
hdfs dfs -copyFromLocal $PATH_TO_FILES$INPUTFILE

hdfs dfs -rm -r -f $OUTPUTFILE
yarn jar "$ABSOLUTE_PATH_TO_REPOSITORY"/peacok_mapred_sort/target/sort-mapreduce-0.0.1.jar "$HDFS_PATH$INPUTFILE" "$HDFS_PATH$OUTPUTFILE" $NB_LINES $NB_REDUCERS

rm -r -f $PATH_TO_FILES$OUTPUTFILE
hdfs dfs -getmerge $OUTPUTFILE $PATH_TO_FILES$OUTPUTFILE

hdfs dfs -rm -r -f $INPUTFILE
hdfs dfs -rm -r -f $OUTPUTFILE
