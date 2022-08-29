#!/bin/bash

#DESCRIPTION:
#bash visualize-UCSC.sh [inputFile] [outputFile] [QvalThresh]
#        
#        [inputFile]                Input Fit-Hi-C file to visualize 
#        [outputFile]               Output file for UCSC to visualize 
#        [QvalThresh]               Q-value threshold to filter Fit-Hi-C interactions at 



INPUT=$1
OUTPUT=$2
QVALTHRESH=$3

echo "track type=interact name=\"Your_Fit-Hi-C_Interactions\" description=\"Fit-Hi-C_Interactions\" interactDirectional=true useScore=on maxHeightPixels=50:100:200 visibility=full" > $OUTPUT
echo "#chrom  chromStart  chromEnd  name  score  value  exp  color  sourceChrom  sourceStart  sourceEnd  sourceName  sourceStrand  targetChrom  targetStart  targetEnd  targetName  targetStrand" >> $OUTPUT
zcat $INPUT | awk -v q="$QVALTHRESH" '{if($7<q){print $0}}' | awk '{print $1, ($2-1), ($4+1), NR, int(-log($7)/log(10)), -log($7)/log(10), "EXP", "0", $1, ($2-1), ($2+1), "SOURCE_NAME", ".", $3, ($4-1), ($4+1), "TARGET_NAME", "+"}' >> $OUTPUT
#max=$(cat $OUTPUT.tmp | awk '{if(NR>2){print $0}}' | sort -nrk5,5 | head -1 | cut -d" " -f5)
#scaling_factor=$(echo "1000/$max" | bc -l)
#cat $OUTPUT.tmp | awk -v sf="$scaling_factor" '{if(NR<2){print $0}else{print $1, $2, $3, $4, int($5*sf), $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18}}' > $OUTPUT
#rm $OUTPUT.tmp
