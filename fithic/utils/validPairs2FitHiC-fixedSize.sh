#!/bin/bash 
set -o nounset
set -o pipefail
set -o errexit 

#DESCRIPTION:
#bash validPairs2FitHiC-fixedSize.sh [resolution] [libraryName] [validPairsFile]        
#[resolution]         The resolution of the dataset being studied
#[libraryName]        The prefix of the file generated
#[validPairsFile]     A textfile containing the validPairs


w=$1
libName=$2
validPairsFile=$3
echo "The resolution given is $w" 
echo "Library name is $libName"
outdir=$4
lowDistThres=$(echo "$w*2" | bc)
echo "The lower distance threshold is $lowDistThres"

################################
## Outline of each awk line ##
    #eliminate random chrs etc -- MAY change from organism to organism
    #eliminate very close intrachr counts -- MAY change from lib to lib, and org to org
    #bin using the window size 
    #sort-uniq to count 
    #print
    #compress
################################

zcat -f $validPairsFile | 
awk '{if(length($2)<=5 && length($5)<=5) {print $0}}' | grep -v chrM | 
awk -v t=$lowDistThres '{if($2!=$5||($2==$5 && (sqrt(($3-$6)^2>t)))) {print $0}}' | 
awk -v r=$w '{print $2, int($3/r)*r, $5, int($6/r)*r}' | 
awk 'BEGIN{OFS="\t"} {if($1<$3 || ($1==$3 && $2<=$4)) {print $1, $2, $3, $4} else {print $3, $4, $1, $2}}' | 
sort | uniq -c |  sed -r 's/^( *[^ ]+) +/\1\t/' | 
awk -v r=$w 'BEGIN{OFS="\t"; FS="\t"} {print $2, $3+r/2, $4, $5+r/2, $1}' |
gzip > $outdir/$libName\_fithic.contactCounts.gz
exit 
