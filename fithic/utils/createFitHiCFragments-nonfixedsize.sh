#!/bin/bash

#DESCRIPTION:
#bash createFitHiCFragments-nonfixedsize.sh [outputFile] [RE] [fastaReferenceGenome]
#        
#        [outputFile]               A desired output file path
#        [RE]                       Either the name of the restriction enzyme used, or the cutting position using “^”. For example, A^AGCTT for HindIII.
#        [fastaReferenceGenome]     A reference genome in fasta format


outfile=$1
cutsite=$2
fasta=$3

wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/bin/utils/digest_genome.py
python2 ./digest_genome.py $fasta -r $cutsite -o ./bedfile.txt
rm ./digest_genome.py
currChr=$(head -n 1 bedfile.txt|cut -f1)
ctr=1
while read p; do
    line=( $p )
    if [ $currChr = ${line[0]} ]; then
        echo -e "chr$ctr\t0\t${line[2]}\t1\t1" >> $outfile
    else
        echo "Working on $currChr"
        ctr=$((ctr+1))
        currChr=${line[0]}
        echo -e "chr$ctr\t0\t${line[2]}\t1\t1" >> $outfile
    fi
done <bedfile.txt
echo "Working on $currChr"
rm bedfile.txt
gzip $outfile
