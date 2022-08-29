#!/bin/bash

HIC_DUMPED=$1
CHR1=$2
CHR2=$3
OUT=$4
RES=$5
if [[ "$#" -gt 5 ]]; then
	NORM=$6
else
	NORM="KR"
fi
if [[ "$#" -gt 6 ]]; then
	DATATYPE=$7
else
	DATATYPE="observed"
fi

currscriptdir=`dirname $0`
echo "currscriptdir : "$currscriptdir
cd $currscriptdir

~/packages/Anaconda3/anaconda3/bin/python3 createFitHiCContacts-hic.py --HiCFile ${HIC_DUMPED} --CHR1 ${CHR1} --CHR2 ${CHR2} --resolution ${RES} --datatype ${DATATYPE} --Norm ${NORM} --outFile ${OUT}

# cat $HIC_DUMPED | awk -v chr1=$CHR1 -v chr2=$CHR2 '{ printf "%s\t%s\t%s\t%s\t%s\n", chr1,$1,chr2,$2,$3}' | gzip > $OUT

