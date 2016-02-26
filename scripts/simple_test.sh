#!/bin/bash

################################################################################
# we download a sample bank from EBI
################################################################################
# if wget is not installed, you may use "curl -O ..."
DATA_SAMPLE="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR387/SRR387476/SRR387476.fastq.gz"
WGET_PATH=`which wget`
echo ">>> Retrieving data sample: ${DATA_SAMPLE}"
if [ ! -z "$WGET_PATH" ] ; then
  echo "    using '$WGET_PATH'..."
  #wget --quiet: to discard progress monitor
  wget ${DATA_SAMPLE}
else
   CURL_PATH=`which curl`
  if [ ! -z "$CURL_PATH" ] ; then
    echo "    using '$CURL_PATH'..."
    #curl --silent: to discard progress monitor
    curl -O ${DATA_SAMPLE}
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

echo ">>> step 2/3: comparing original file and leon processed one..."
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
#rm -f  SRR387476.*

exit $RETVAL