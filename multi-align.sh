# Joseph Zambreno
# Addtions by Mihir Awatramani
# multi-align.sh - Script that runs the MEMOCODE 2012 challenge by partitioning 
# the sequence data evenly across the specified list of machines.


#!/bin/bash


# Local hostname, date, etc.
MYNAME=`hostname`
MYDATE=`date +%F`
MYTIME=`date +%c`
MYUSER=`whoami`

# Output directory. For now, just use the username information
OUTDIR="/var/tmp/"$MYUSER"-memocode/"
OUT="$OUTDIR'/'output.txt"

# We are specific with the command-line arguments
NUMARGS=$#
if [ $NUMARGS -lt 3 ]; then
	echo "Usage: $0 file.lst genome.bin sequences.bin start_seq end_seq. The start and end sequence is not specified they would default to start and end of the file."
	exit;
fi

# Machines to run on, provided as a list in a file
MACHINES=`cat $1`
# Human genome file
GENOME=$2
# Subsequence file
SEQUENCES=$3

if [ $NUMARGS -ge 4 ];then
# start sequence value
 START_SEQ=$4
if [ $NUMARGS -ge 5 ];then
# end sequence value
 END_SEQ=$5
else
 FILESIZE=$(stat -c%s "$3")
 END_SEQ=$[ $FILESIZE/32 ]
fi
else
START_SEQ=0
fi 


# Run the aligner application. We do not cat all the results. So we wait for the application to finish before moving to the next machine
echo "Running the aligner application"

for i in $MACHINES
do

  echo "Running ./hash-align $GENOME $SEQUENCES 100 $START $END on $i..."
 
  ssh -o StrictHostKeyChecking=no $i "export PATH=/usr/local/cuda/bin:$PATH;export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH;cd $OUTDIR;./hash-align $GENOME $SEQUENCES 100 $i $START_SEQ $END_SEQ > output3_on_$i.txt &"

done




