#!/bin/bash

# One can define SILENT_MODE=true to enable progress monitor while loading data
# from EBI.

################################################################################
# we download a sample bank from EBI
################################################################################
# if wget is not installed, you may use "curl -O ..."
rm -f  SRR387476.*
DATA_SAMPLE="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR387/SRR387476/SRR387476.fastq.gz"
WGET_PATH=`which wget`
echo ">>> Retrieving data sample: ${DATA_SAMPLE}"
if [ ! -z "$WGET_PATH" ] ; then
  echo "    using '$WGET_PATH'..."
  if [ $SILENT_MODE=="true"  ] ; then
    wget --quiet ${DATA_SAMPLE}
  else
    wget ${DATA_SAMPLE}
  fi
else
   CURL_PATH=`which curl`
  if [ ! -z "$CURL_PATH" ] ; then
    echo "    using '$CURL_PATH'..."
    if [ $SILENT_MODE=="true"  ] ; then
      curl --silent -O ${DATA_SAMPLE}
    else
      curl -O ${DATA_SAMPLE}
    fi
  else
    echo "    /!\ error: unable to find 'wget' or 'curl'"
    exit 1
  fi
fi

echo
echo ">>> Gunzip data sample..."
gunzip SRR387476.fastq.gz

################################################################################
# we launch leon
################################################################################
echo
echo ">>> Start 'leon': "

if [ -e ../build/bin/leon ]; then
  # when working in devel mode
  LEON_CMD="../build/bin/leon"
else
  # when working in prod mode
  LEON_CMD="./bin/leon"
fi

echo ">>> step 1/3: running compressing process..."
${LEON_CMD} -c -lossless -file SRR387476.fastq
echo
echo ">>> step 2/3: running uncompressing process..."
${LEON_CMD} -d -file SRR387476.fastq.leon

################################################################################
# we check the result
################################################################################

echo ">>> step 3/3: comparing original file and leon processed one..."
diff ./SRR387476.fastq ./SRR387476.fastq.d
if [ $? -eq 0 ]; then
   echo "*** test OK ***  "
   RETVAL=0
else
   echo "/!\/!\  test KO  /!\/!\  "
   RETVAL=1
fi

################################################################################
# clean up
################################################################################
rm -f  SRR387476.*

exit $RETVAL